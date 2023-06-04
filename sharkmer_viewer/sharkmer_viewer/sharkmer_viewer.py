# Ingests a tsv file of kmer histograms, where the first column is the count
# of kmers. The following columns are the kmer counts for each sample.
# Outputs an animated histogram of the kmer counts for each sample.

import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import scipy
import plotly.graph_objects as go
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler

peak_threshold = 1000

def get_limits(df_histo):
    # Calculate the limits of the plot based on the characteristics of the peak, if there is one
    # Get the last column of the dataframe as a numpy array
    y = df_histo.iloc[:, -1]
    y = np.array(y)

    x_limit = 100
    y_limit = 50

    y_max = None
    # Loop through the columns of the dataframe and find the tallest peak, use its height to scale y
    for i in range(len(df_histo.columns)):
        y = df_histo.iloc[:, i]
        y = np.array(y)

        # Find the peaks
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        peaks, properties = scipy.signal.find_peaks(y, distance=2, threshold=peak_threshold )

        x_of_tallest_peak = None
        for peak in peaks:
            # Ignore peaks that are too far to the left
            if peak > 5:
                if y_max is None:
                    y_max = y[peak]
                    x_of_tallest_peak = peak
                elif y[peak] > y_max:
                    y_max = y[peak]
                    x_of_tallest_peak = peak
        
        # If last column
        if i == len(df_histo.columns) - 1:
            if x_of_tallest_peak is not None:
                x_limit = x_of_tallest_peak * 3

    if y_max is not None:
        y_limit = y_max * 1.2
    
    return x_limit, y_limit

def get_tallest_peaks(y):
    # Get a vector of the peaks in descending order of height
    peaks, _ = scipy.signal.find_peaks(y, threshold=peak_threshold)
    # Sort the peaks by height
    peaks = sorted(peaks, key=lambda x: y[x], reverse=True)
    return peaks

def reindex(df):
    # Create a new column that is the index of the feature, where the index is the order of the feature
    # in the last sample where both features are present

    # Get the last sample with all feature indexes
    all_raw_indexes = df['raw_index'].unique()

    if len(all_raw_indexes) < 2:
        df['index'] = df['raw_index'].astype(int)
        return df

    all_samples = df['sample'].unique()

    # Get the last sample with all raw_indexes
    last_sample_with_all_indexes = None
    for sample in all_samples:
        df_sub = df[df['sample'] == sample]
        raw_indexes = df_sub['raw_index'].unique()
        if len(raw_indexes) == len(all_raw_indexes):
            if last_sample_with_all_indexes is None:
                last_sample_with_all_indexes = sample
            elif sample > last_sample_with_all_indexes:
                last_sample_with_all_indexes = sample
    
    if last_sample_with_all_indexes is None:
        # Can't order the features because there is no sample with all features
        df['index'] = df['raw_index'].astype(int)
        return df
    
    # Create a new index numbering scheme, where the index is the order of the coverage in the last sample
    df_sub = df[df['sample'] == last_sample_with_all_indexes]
    df_sub = df_sub.sort_values(by=['coverage'])
    df_sub['index'] = np.arange(len(df_sub))

    # update the index column in the original dataframe with the new index
    df = df.merge(df_sub[['raw_index', 'index']], how='left', on='raw_index')
    df['index'] = df['index'].astype(int)

    return df


def create_report(in_histo_name, in_stats_name, out_name, run_name, genome_size):
    # Read in the histogram data
    df_histo = pd.read_csv(in_histo_name, sep="\t", header=None)

    # Get the counts, then remove the column
    x = df_histo.iloc[:, 0]
    x = np.array(x)
    df_histo = df_histo.drop(df_histo.columns[0], axis=1)

    # Truncate the data
    df_histo = df_histo.iloc[:100]

    # Parse the stats file from sharkmer
    stats_dict = {}
    with open(in_stats_name, "r") as f:
        for line in f:
            line = line.strip()
            line = line.split("\t")
            stats_dict[line[0]] = line[1]

    # Create a vector with the cumulative bases read for each sample
    n_samples = len(df_histo.columns)
    n_bases_read = int(stats_dict["n_bases_read"])
    n_bases_per_sample = n_bases_read // n_samples
    cumulative_bases_read = np.arange(n_bases_per_sample, n_bases_read + 1, n_bases_per_sample)
    cumulative_coverage = cumulative_bases_read / 1000000 / genome_size

    # Get the number of peaks in the final column
    y = df_histo.iloc[:, -1]
    y = np.array(y)
    peaks = get_tallest_peaks(y)
    n_peaks = len(peaks)

    if n_peaks == 0:
        print("No peaks found")
        return 0
    elif n_peaks > 2:
        print("More than two peaks found, for now we only support diploid genomes")
        return 0

    # Create a new data frame of peaks. 
    df_peaks = pd.DataFrame(columns=["sample", "coverage", "frequency"])
    for i in range(len(df_histo.columns)):
        y = df_histo.iloc[:, i]
        y = np.array(y)
        peaks, _ = scipy.signal.find_peaks(y, threshold=peak_threshold)
        for peak in peaks:
            df_new_row = pd.DataFrame({"sample": [i], "coverage": [peak], "frequency": [y[peak]]})
            df_peaks = pd.concat([df_peaks, df_new_row], ignore_index=True)

    # Use spectral clustering on peak_index and peak_height to cluster the peaks
    # https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html
    scaler = StandardScaler()
    df_peaks_scaled = scaler.fit_transform(df_peaks[['coverage', 'frequency']])
    clustering = SpectralClustering(n_clusters=n_peaks, assign_labels="discretize", random_state=0).fit(df_peaks_scaled)
    df_peaks['feature'] = "peak"
    df_peaks['raw_index'] = clustering.labels_

    df_peaks = reindex(df_peaks)

    # Create a new data frame of valleys.
    df_valleys = pd.DataFrame(columns=["sample", "coverage", "frequency"])
    for i in range(len(df_histo.columns)):
        y = df_histo.iloc[:, i]
        y = np.array(y)
        valleys, _ = scipy.signal.find_peaks(-y, threshold=peak_threshold)
        for valley in valleys:
            df_new_row = pd.DataFrame({"sample": [i], "coverage": [valley], "frequency": [y[valley]]})
            df_valleys = pd.concat([df_valleys, df_new_row], ignore_index=True)
    
    # Use spectral clustering on valley_index and valley_height to cluster the valleys
    # https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html
    scaler = StandardScaler()
    df_valleys_scaled = scaler.fit_transform(df_valleys[['coverage', 'frequency']])
    clustering = SpectralClustering(n_clusters=n_peaks, assign_labels="discretize", random_state=0).fit(df_valleys_scaled)
    df_valleys['feature'] = "valley"
    df_valleys['raw_index'] = clustering.labels_

    df_valleys = reindex(df_valleys)

    # Combine the peaks and valleys into a single dataframe
    df_features = pd.concat([df_peaks, df_valleys], ignore_index=True)
    df_features['name'] = df_features['feature'] + "_" + df_features['index'].astype(str)

    # Add a column with feature_symbols based on feature column for the plot, where peaks are triangles and valleys are circles
    df_features['feature_symbol'] = df_features['feature']
    df_features.loc[df_features['feature'] == 'peak', 'feature_symbol'] = 'triangle-up'
    df_features.loc[df_features['feature'] == 'valley', 'feature_symbol'] = 'circle'

    # Plot the features and lines
    # https://plotly.com/python/line-and-scatter/
    duration = 50 # milliseconds
    unique_names = df_features['name'].unique()

    # First plot is created with the scatter points
    fig_features = go.Figure()

    for name in unique_names:
        df_sub = df_features[df_features['name'] == name]
        fig_features.add_trace(go.Scatter(
            x=df_sub['coverage'], 
            y=df_sub['frequency'],
            mode='markers', 
            marker_symbol=df_sub['feature_symbol'], 
            name=name,
            text=df_sub['sample'],  
            hovertemplate = 'Coverage: %{x}<br>Frequency: %{y}<br>Sample: %{text}'  
        ))

    # Second plot is created with line traces
    histo_color = 'rgba(86, 180, 233, 0.5)'
    fig_histo = go.Figure(
        data=[go.Scatter(x=x, y=[0]*len(x), mode='lines')],
        layout=go.Layout(
            xaxis=dict(range=[0, get_limits(df_histo)[0]], autorange=False, title='Coverage'),
            yaxis=dict(range=[0, get_limits(df_histo)[1]], autorange=False, title='Frequency'),
            updatemenus=[dict(type="buttons",
                            buttons=[dict(label="Play",
                                            method="animate",
                                            args=[None, {"frame": {"duration": duration, "redraw": True}, 
                                                        "fromcurrent": True, 
                                                        "transition": {"duration": duration, "easing": "cubic-in-out"}}])])],
            annotations=[dict(x=1, y=1, xref='paper', yref='paper', text=run_name, showarrow=False, font=dict(size=20))]),
        frames=[go.Frame(
            data=[go.Scatter(
                x=x,
                y=df_histo.iloc[:, i],
                mode='lines',
                fill='tozeroy',
                fillcolor=histo_color,
                marker=dict(color=histo_color, size=10))])
            for i in range(len(df_histo.columns))]
    )

    # Adding traces from fig_features to fig_histo
    for trace in fig_features['data']:
        fig_histo.add_trace(trace)

    # Now fig_histo contains both the line traces and scatter points
    fig_histo.show()
    fig_histo.write_html(out_name + ".html")

    # If there are two peaks, calculate the genome size with the manual method
    # https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
    if n_peaks == 2:
        df_estimates = pd.DataFrame(columns=[
            "sample", 
            "cumulative_bases_read", 
            "first_valley", 
            "heterozygous_peak", 
            "homozygous_peak", 
            "genome_size_heterozygous", 
            "genome_size_homozygous"])

        # Populate df_estimates with sample, cumulative_bases_read, and NaN for the rest
        for i in range(len(df_histo.columns)):
            df_new_row = pd.DataFrame({"sample": [i+1], "cumulative_bases_read": [cumulative_bases_read[i]], "first_valley": [np.nan], "heterozygous_peak": [np.nan], "homozygous_peak": [np.nan]})
            df_estimates = pd.concat([df_estimates, df_new_row], ignore_index=True)
        
        # Populate df_estimates['first_valley'] with the valleys with index 0
        for i in range(len(df_histo.columns)):
            df_sub = df_features[
                (df_features['sample'] == i) & 
                (df_features['feature'] == 'valley') & 
                (df_features['index'] == 0)]
            for sample, coverage in zip(df_sub['sample'], df_sub['coverage']):
                df_estimates.loc[df_estimates['sample'] == sample, 'first_valley'] = coverage

        # Populate df_estimates['heterozygous_peak'] with the peaks with index 0
        for i in range(len(df_histo.columns)):
            df_sub = df_features[
                (df_features['sample'] == i) & 
                (df_features['feature'] == 'peak') & 
                (df_features['index'] == 0)]
            for sample, coverage in zip(df_sub['sample'], df_sub['coverage']):
                df_estimates.loc[df_estimates['sample'] == sample, 'heterozygous_peak'] = coverage

        # Populate df_estimates['homozygous_peak'] with the peaks with index 1
        for i in range(len(df_histo.columns)):
            df_sub = df_features[
                (df_features['sample'] == i) & 
                (df_features['feature'] == 'peak') & 
                (df_features['index'] == 1)]
            for sample, coverage in zip(df_sub['sample'], df_sub['coverage']):
                df_estimates.loc[df_estimates['sample'] == sample, 'homozygous_peak'] = coverage

        def integrate_histo_kmers(df_histo, column, end):
            histo = df_histo.iloc[:, column]
            histo = np.array(histo)
            integral = 0
            for i in range(end):
                # Multiply the number of kmers by the coverage, which is the index plus 1
                integral += histo[i] * (i+1)
            return integral

        # Calculate the genome sizes
        for i in range(len(df_histo.columns)):

            if not np.isnan(df_estimates["first_valley"][i]):
                histo = df_histo.iloc[:, i]
                n_kmers_all = integrate_histo_kmers(df_histo, i, len(df_histo))
                n_kmers_error = integrate_histo_kmers(df_histo, i, int(df_estimates['first_valley'][i]))
                n_kmers = n_kmers_all - n_kmers_error

                if not np.isnan(df_estimates["heterozygous_peak"][i]):
                    heterozygous_peak = df_estimates['heterozygous_peak'][i]
                    df_estimates["genome_size_heterozygous"][i] = n_kmers / 2 / heterozygous_peak

                if not np.isnan(df_estimates.iloc[i, 4]):
                    homozygous_peak = df_estimates['homozygous_peak'][i]
                    df_estimates["genome_size_homozygous"][i] = n_kmers / homozygous_peak

    # print all lines of the dataframe
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df_estimates)
    
    # Plot the heterozygous and homozygous genome size estimates for each sample
    fig_genome_size = go.Figure()
    fig_genome_size.add_trace(go.Scatter(
        x=df_estimates['sample'],
        y=df_estimates['genome_size_heterozygous']/1e6,
        mode='lines',
        name='heterozygous genome size estimate'
    ))
    fig_genome_size.add_trace(go.Scatter(
        x=df_estimates['sample'],
        y=df_estimates['genome_size_homozygous']/1e6,
        mode='lines',
        name='homozygous genome size estimate'
    ))
    fig_genome_size.update_layout(
        xaxis_title="Sample",
        yaxis_title="Genome size (Mb)",
        title=run_name
    )
    fig_genome_size.show()
    fig_genome_size.write_html(out_name + "_genome_size.html")

    return 0


def main():
    # https://docs.python.org/3/howto/argparse.html
    parser = argparse.ArgumentParser(description="view sharkmer results")
    parser.add_argument(
        "-d", "--histogram",
        help="input histogram distribution file from sharkmer",
        required=True
    )

    parser.add_argument(
        "-s", "--stats",
        help="input stats file from sharkmer",
        required=True
    )

    parser.add_argument(
        "-n", "--name",
        help="run name used in output files, optional",
        default=""  # Default value if no output file is specified
    )

    parser.add_argument(
        "-o", "--output",
        help="output file base name, optional",
        default=""  # Default value if no output file is specified
    )

    parser.add_argument(
        "-g", "--genome-size",
        help="haploid genome size in megabases, optional",
        type=float,
        default=None  # Default value if no output file is specified
    )  

    args = parser.parse_args()

    in_histo_name = args.histogram
    run_name = args.name
    in_stats_name = args.stats
    genome_size = args.genome_size

    # If no output file is specified, use the histo name without .histo at the end
    out_name = args.output
    if out_name == "":
        out_name = in_histo_name.replace(".histo", "")


    create_report(in_histo_name, in_stats_name, out_name, run_name, genome_size)


if __name__ == "__main__":
    main()
