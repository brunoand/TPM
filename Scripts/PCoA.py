#!/usr/bin/python3

import skbio
import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import warnings
import sys

Input = sys.argv[1]
Meta = sys.argv[2]
Output = sys.argv[3]
fig = sys.argv[4]

warnings.filterwarnings("ignore")
metadata = pd.read_csv(Meta, sep = '\t', index_col = 0)
my_obj =  DistanceMatrix.read(Input, 'lsmat')
PC = pcoa(my_obj)



def plot_PCoA(matrix, ID_column, fig):
    figure = PC.plot(metadata, ID_column, axis_labels=('PC 1 explains {0:.2f}%'.format(PC.proportion_explained[0]*100), 'PC 2 explains {0:.2f}%'.format(PC.proportion_explained[1]*100), 'PC 3 explains {0:.2f}%'.format(PC.proportion_explained[2]*100)), cmap='jet', s=50)
    figure.set_size_inches(12.5, 8.5)
    figure.text(0,0.9, r'Samples colored by {}'.format(ID_column), fontsize=16)
    figure.savefig(Output + 'PCOA_{}.'.format(ID_column) + fig, bbox_inches='tight')
    
for x in metadata.columns:
    plot_PCoA(PC, x, fig)
