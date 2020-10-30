import os
import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wasserstein_distance



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
sns.set(style='dark')
font={'size':10 }# 'family':'serif','style':'italic','weight':'normal','color':'red',
      

plt.figure(figsize=(7, 4.8)) # default 6.4, 4.8
###########IFG section###########
dataset = pd.read_csv('coverage_train_dataset.csv')
#ifg=dataset[dataset["type"]=='IFG']
sns.barplot(x='range',y='coverage', hue='type',data=dataset, hue_order=['RS','FG'])
plt.xlabel('Frequency of Occurrence in GDB13')
plt.ylabel('Coverage (%)')
title='The GDB13 coverage of 1M training dataset'
sns.despine(bottom=True)
plt.title(title, fontsize=14)
plt.legend()
ifg=dataset[dataset['type']=='RS']
for x, y in enumerate(ifg['coverage'].values):
    plt.text(x-0.45, round(y,2), "%.2f"%y, fontdict=font)
#plt.tight_layout()
ifg=dataset[dataset['type']=='FG']
for x, y in enumerate(ifg['coverage'].values):
    plt.text(x, round(y,2), "%.2f"%y, fontdict=font)
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title+'.png'),
    dpi=250
)


