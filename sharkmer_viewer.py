# Ingests a tsv file of kmer histograms, where the first column is the count 
# of kmers. The following columns are the kmer counts for each sample.
# Outputs an animated histogram of the kmer counts for each sample.

import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t', header=None)
df = df.iloc[:100]

fig, ax = plt.subplots()

x = df['x']
lines = []
for column in df.columns[1:]:
    line, = ax.plot(x, df[column])
    lines.append(line)

def animate(i):
    for j, line in enumerate(lines):
        line.set_ydata(df[df.columns[j+1]] * np.sin(i / 100.))  # update the data
    return lines

# Init only required for blitting to give a clean slate.
def init():
    for line in lines:
        line.set_ydata(np.ma.array(x, mask=True))
    return lines

ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
                              interval=25, blit=True)

# Save the animation
# ani.save('animation.mp4', writer='ffmpeg')
ani.save('animation.gif', writer='imagemagick')
# plt.show()