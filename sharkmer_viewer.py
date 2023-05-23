# Ingests a tsv file of kmer histograms, where the first column is the count 
# of kmers. The following columns are the kmer counts for each sample.
# Outputs an animated histogram of the kmer counts for each sample.

import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

histo = pd.read_csv(sys.argv[1], sep='\t', header=None)

