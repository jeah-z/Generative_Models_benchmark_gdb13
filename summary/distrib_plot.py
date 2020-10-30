import os
import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wasserstein_distance

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2.-0.2, 1.03*height, '%s' % float(height))

def get_parser():
    parser = argparse.ArgumentParser(
        "Prepares distribution plots for weight, logP, SA, and QED\n"
    )
    parser.add_argument(
        '--test', type=str, default='test.csv', help='Path to the test set'
    )
    parser.add_argument(
        '--config', '-c', type=str, default='config.csv',
        help='Path to the config csv with `name` and `path` columns. '
             '`name` is a model name, and '
             '`path` is a path to generated samples`'
    )
    parser.add_argument(
        '--n_jobs', type=int, default=1,
        help='number of processes to use'
    )
    parser.add_argument(
        '--img_folder', type=str, default='images/',
        help='Store images in this folder'
    )
    return parser
parser = get_parser()
config, unknown = parser.parse_known_args()
sns.set(style='ticks')

plt.figure(figsize=(7.4, 4.8)) # default 6.4, 4.8
dataset = pd.read_csv('gdb_distrib.csv')

######## GDB13 distribution ########
plt.clf()
ifg=dataset
g=sns.barplot(x='range',y='percentage', hue='Type',data=ifg,hue_order=['RS','FG'])
plt.xlabel('Frequency of Occurrence in GDB-13')
plt.yscale('log')
plt.ylabel('Percentage (%)')
title='Distribution of ring systems and functional groups in GDB-13'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0 )
plt.title(title, fontsize=14)
plt.legend()
# ifg=ifg[ifg['Type']=='RS']
# for index, row in ifg.iterrows():
#     g.text(index*5, row.percentage, row.percentage, color='black', ha="center")
# ifg=dataset[dataset['Type']=='RS']
# for x, y in enumerate(ifg['percentage'].values):
#     plt.text(x-0.45, round(y,2), "%.2f"%y)
# #plt.tight_layout()
# ifg=dataset[dataset['Type']=='FG']
# for x, y in enumerate(ifg['percentage'].values):
#     plt.text(x, round(y,2), "%.2f"%y)
#plt.show()
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.png'),
    dpi=250
)