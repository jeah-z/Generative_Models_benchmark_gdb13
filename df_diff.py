from rdkit import Chem
from collections import namedtuple
import rdkit
import numpy as np
import pandas as pd
from functools import partial
from multiprocessing import Pool
import pandas as pd
from tqdm.auto import tqdm
import pathlib
import argparse


def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input_file', '-i', type=str,
        default='./',
        help='path to sampled files'
    )
    parser.add_argument(
        '--target_file', '-t', type=str,
        default='./SampledDatasetCano.csv',
        help='path to save files'
    )

    return parser


def pd_diff(input, target, config):
    same_flag = input['SMILES'].isin(target['SMILES'])
    diff_flag = [not f for f in same_flag]
    same = input[same_flag]
    diff = input[diff_flag]
    same.to_csv(config.input_file+'_same.csv', index=None)
    diff.to_csv(config.input_file+'_diff.csv', index=None)


def main(config):
    file = config.input_file+'.csv'
    target = config.target_file+'.csv'
    file_pd = pd.read_csv(file, names=['SMILES', 'model'], low_memory=False)
    print('%s was loaded!\n' % (file))
    target_pd = pd.read_csv(
        target, names=['SMILES', 'model'], low_memory=False)
    print('%s was loaded!\n' % (target))
    pd_diff(file_pd, target_pd, config)


if __name__ == "__main__":
    parser = get_parser()
    config, unknown = parser.parse_known_args()
    main(config)
