#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
plt.switch_backend('agg') 
import sys



def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j, regex=True)
    return text


def parse(Dataframe, Output, Percentage, Fig):
    taxa_all = Dataframe.reset_index()
    taxa_all = taxa_all.drop(['Score','OTU'], axis=1)
    taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= 's__.{1,}', value='', regex=True)
    taxa_all['Taxonomy'] = taxa_all['Taxonomy'].replace(to_replace= '(.__|\s.{1,}|\[|\]|\-.{1,})', value='', regex=True)
    taxa_all2 = pd.DataFrame(taxa_all['Taxonomy'].str.split(';').values.tolist(), columns = ['Domain', 'Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'], index = taxa_all.index)

    taxa_all2 = taxa_all2.drop('Species', axis=1)
    taxa_all2 = taxa_all2.merge(taxa_all, left_index=True, right_index = True).drop('Taxonomy', axis=1)
    for x in ['Phylum', 'Order', 'Class', 'Family', 'Genus']:
        rel_abundance(taxa_all2, x, Output, Percentage, Fig)


Input_OTU = sys.argv[1]
Taxonomy = sys.argv[2]
Output = sys.argv[3]

OTU_table = pd.read_csv(Input_OTU, sep = '\t')
header_taxonomy = ['OTU', 'Taxonomy', 'Score', 'Size']
Taxonomy = pd.read_csv(Taxonomy, sep = '\t', names = header_taxonomy)
Pattern = {"D_0":"d", "D_1":"p", "D_2":"c", "D_3":"o", "D_4":"f", "D_5":"g", "D_6":"s" }

Taxonomy['Taxonomy']= replace_all(Taxonomy['Taxonomy'], Pattern)
Taxonomy = Taxonomy.merge(OTU_table, right_on='OTU', left_on = 'OTU').drop(['Size'], axis=1).set_index('OTU')
Taxonomy.to_csv(Output, sep = '\t')

