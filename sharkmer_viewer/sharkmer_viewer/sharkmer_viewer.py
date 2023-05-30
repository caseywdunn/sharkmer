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

def get_limits(df):
    # Calculate the limits of the plot based on the characteristics of the peak, if there is one
    # Get the last column of the dataframe as a numpy array
    y = df.iloc[:, -1]
    y = np.array(y)

    x_limit = 100
    y_limit = 50

    y_max = None
    # Loop through the columns of the dataframe and find the tallest peak, use its height to scale y
    for i in range(len(df.columns)):
        y = df.iloc[:, i]
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
        if i == len(df.columns) - 1:
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


def create_report(in_file_name):
    df = pd.read_csv(in_file_name, sep="\t", header=None)

    # Truncate the data
    df = df.iloc[:100]

    # Get the counts, then remove the column
    x = df.iloc[:, 0]
    x = np.array(x)
    df = df.drop(df.columns[0], axis=1)

    # Create the plot
    fig = go.Figure(
        data=[go.Scatter(x=x, y=[0]*len(x), mode='lines')],
        layout=go.Layout(
            xaxis=dict(range=[0, get_limits(df)[0]], autorange=False),
            yaxis=dict(range=[0, get_limits(df)[1]], autorange=False),
            updatemenus=[dict(type="buttons",
                              buttons=[dict(label="Play",
                                            method="animate",
                                            args=[None, {"frame": {"duration": 200, "redraw": True}, 
                                                         "fromcurrent": True, 
                                                         "transition": {"duration": 200, "easing": "cubic-in-out"}}])])]),
        frames=[go.Frame(
            data=[go.Scatter(
                x=x,
                y=df.iloc[:, i],
                mode='lines',
                marker=dict(color="red", size=10))])
            for i in range(len(df.columns))]
    )

    fig.show()
    fig.write_html(in_file_name + ".html")

    return 0


def main():
    # https://docs.python.org/3/howto/argparse.html
    parser = argparse.ArgumentParser(description="view sharkmer results")
    parser.add_argument(
        "input",
        help="input file, histogram from sharkmer",
    )

    args = parser.parse_args()

    in_file_name = args.input

    create_report(in_file_name)


if __name__ == "__main__":
    main()
