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
sns.set(style='ticks')

plt.figure(figsize=(7, 4.8)) # default 6.4, 4.8
###########Training speed###########
dataset = pd.read_csv('train_speed.csv')
ifg=dataset
sns.barplot(x='Models',y='Train_speed',data=ifg, order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT'])
plt.xlabel('Models')
plt.ylabel('Training speed (min/epoch)')
title='Training speed'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0)
plt.title(title, fontsize=14)
#plt.legend(loc='best',title=title)
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title+'.png'),
    dpi=250
)
###########Sampling speed###########
plt.clf()
dataset = pd.read_csv('train_speed.csv')
ifg=dataset
sns.barplot(x='Models',y='Sample_speed',data=ifg, order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT'])
plt.xlabel('Models')
plt.ylabel('Sampling speed (k compounds/s)')
title = 'Sampling speed'
plt.tick_params(axis=u'x', direction=u'out',length=0 )
sns.despine(bottom=True)
plt.title(title, fontsize=14)
#plt.legend(loc='best',title=title)
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title+'.png'),
    dpi=250
)