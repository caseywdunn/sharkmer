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
    df_histo = pd.read_csv(in_histo_name, sep="\t", header=None)
    # df_stats = pd.read_csv(in_stats_name, sep="\t", header=None)
    # read the tab delimeted stats file and load into a dictionary, where the first column is the key
    # and the second column is the value
    stats_dict = {}
    with open(in_stats_name, "r") as f:
        for line in f:
            line = line.strip()
            line = line.split("\t")
            stats_dict[line[0]] = line[1]

    # Truncate the data
    df_histo = df_histo.iloc[:100]

    # Get the counts, then remove the column
    x = df_histo.iloc[:, 0]
    x = np.array(x)
    df_histo = df_histo.drop(df_histo.columns[0], axis=1)

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
    fig1 = go.Figure()

    for name in unique_names:
        df_sub = df_features[df_features['name'] == name]
        fig1.add_trace(go.Scatter(
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
    fig2 = go.Figure(
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

    # Adding traces from fig1 to fig2
    for trace in fig1['data']:
        fig2.add_trace(trace)

    # Now fig2 contains both the line traces and scatter points
    fig2.show()
    fig2.write_html(out_name + ".html")



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
        help="output file, optional",
        default=""  # Default value if no output file is specified
    )

    parser.add_argument(
        "-o", "--output",
        help="output base name, optional",
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
