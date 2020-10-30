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

plt.figure(figsize=(8, 4.8)) # default 6.4, 4.8
###########IFG section###########
dataset = pd.read_csv('coverage.csv')
ifg=dataset[dataset["type"]=='IFG']
sns.barplot(x='range',y='coverage', hue='Model',data=ifg, hue_order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT','Training_set'])
plt.xlabel('Frequency of Occurrence in GDB-13')
plt.ylabel('Percentage (%)')
title='Coverage of Functional Groups'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0)
plt.title(title, fontsize=14)
plt.legend()
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title+'.png'),
    dpi=250
)
########RS section########
plt.clf()
ifg=dataset[dataset["type"]=='RS']
sns.barplot(x='range',y='coverage', hue='Model',data=ifg,hue_order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT','Training_set'])
plt.xlabel('Frequency of Occurrence in GDB-13')
plt.ylabel('Percentage (%)')
title='Coverage of Ring Systems'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0)
plt.title(title, fontsize=14)
plt.legend(bbox_to_anchor=(0.80, 0.55))
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title+'.png'),
    dpi=250
)

######## IFG generated section########
plt.clf()
ifg=dataset[dataset["type"]=='IFG']
sns.barplot(x='range',y='percentage_model', hue='Model',data=ifg,hue_order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT'])
plt.xlabel('Frequency of Occurrence in GDB-13')
plt.ylabel('Percentage (%)')
plt.yscale('log')
title='Distribution of Functional Groups Shared with GDB-13'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0)
plt.title(title, fontsize=14)
plt.legend()
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.png'),
    dpi=250
)

######## RS generated section########
plt.clf()
ifg=dataset[dataset["type"]=='RS']
sns.barplot(x='range',y='percentage_model', hue='Model',data=ifg,hue_order=['REINVENT','CharRNN','AAE','ORGAN','LatentGAN','VAE','GraphINVENT'])
plt.xlabel('Frequency of Occurrence in GDB-13')
plt.ylabel('Percentage (%)')
plt.yscale('log')
title='Distribution of Ring Systems Shared with GDB-13'
sns.despine(bottom=True)
plt.tick_params(axis=u'x', direction=u'out',length=0)
plt.title(title, fontsize=14)
plt.legend()
#plt.tight_layout()
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.pdf')
)
plt.savefig(
    os.path.join(config.img_folder, title.replace('/','_')+'.png'),
    dpi=250
)