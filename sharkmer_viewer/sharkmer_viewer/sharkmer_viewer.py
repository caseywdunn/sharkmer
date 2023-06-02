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
        peaks, properties = scipy.signal.find_peaks(y, height=0, distance=2, threshold=10)

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
    peaks, _ = scipy.signal.find_peaks(y, height=0)
    # Sort the peaks by height
    peaks = sorted(peaks, key=lambda x: y[x], reverse=True)
    return peaks


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

    # Create a new data frame of peaks. The columns are:
    # 1. The sample number (ie the column in df_histo)
    # 2. The peak index
    # 3. The peak height

    df_peaks = pd.DataFrame(columns=["sample", "peak_index", "peak_height"])
    for i in range(len(df_histo.columns)):
        y = df_histo.iloc[:, i]
        y = np.array(y)
        peaks, _ = scipy.signal.find_peaks(y, height=0)
        for peak in peaks:
            df_new_row = pd.DataFrame({"sample": [i], "peak_index": [peak], "peak_height": [y[peak]]})
            df_peaks = pd.concat([df_peaks, df_new_row], ignore_index=True)
    
    # Plot the peaks, where the x axis is the peak index the y axis is the peak height, and the color is the sample
    # https://plotly.com/python/line-and-scatter/
    fig = go.Figure()
    for i in range(len(df_histo.columns)):
        df_peaks_sample = df_peaks[df_peaks["sample"] == i]
        fig.add_trace(go.Scatter(x=df_peaks_sample["peak_index"], y=df_peaks_sample["peak_height"], mode='markers', name="sample " + str(i)))

    fig.update_layout(
        title="Peaks",
        xaxis_title="Peak index",
        yaxis_title="Peak height",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="#7f7f7f"
        )
    )

    fig.show()
    fig.write_html(out_name + "_peaks.html")


    # Create the plot
    duration = 100
    fig = go.Figure(
        data=[go.Scatter(x=x, y=[0]*len(x), mode='lines')],
        layout=go.Layout(
            xaxis=dict(range=[0, get_limits(df_histo)[0]], autorange=False),
            yaxis=dict(range=[0, get_limits(df_histo)[1]], autorange=False),
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
                marker=dict(color="red", size=10))])
            for i in range(len(df_histo.columns))]
    )

    fig.show()
    fig.write_html(out_name + ".html")

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
