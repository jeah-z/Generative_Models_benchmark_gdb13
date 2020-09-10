import rdkit
from rdkit import Chem
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
        '--input', '-i',
        type=str, default='',
        help='the complete name of the input file'
    )
    parser.add_argument(
        '--epoch', '-e',
        type=int, default=100,
        help='epoch that SMILES was sampled'
    )
    parser.add_argument(
        '--n_jobs', '-n',type=int, default=40,
        help='number of processes to use'
    )

    parser.add_argument(
        '--path', type=str,
        default='./',
        help='path to sampled files'
    )
    parser.add_argument(
        '--output', '-o', type=str,
        default='./SampledDatasetCano.csv',
        help='path to save files'
    )
    parser.add_argument(
        '--sample_id', '-si', type=int,
        default='0',
        help='Index of sampling'
    )
    return parser


def file_collection(model, epoch):
    file_set = []
    for i in range(8):
        file_set.append('%s_model_s_%d_%d.csv' % (model, epoch, i))
    return file_set


def open_file(file_set, path):
    dataset = []
    file_path = pathlib.Path(path, file_set)
    with open(file_path, 'r') as f:         # gzip.open(path) as smi:
        next(f)
        lines = f.readlines()
        for line in lines:
            try:
                line = line.replace(',',' ').split()
                dataset.append(line[0])
            except:
                print(str(line) + " is blank line?")
    return  dataset


def atom_count(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        atoms = mol.GetAtoms()
        natoms = len(atoms)
        #smi_cano = Chem.MolToSmiles(mol, isomericSmiles=False)
        return natoms
    except:
        print(smi + "was not valid SMILES\n")
        return None


def counts(n_jobs, output, file_set):
    global dataset
    with Pool(n_jobs) as pool:
        atom_count_p = partial(atom_count)
        dataset_cano = [x for x in tqdm(
            pool.imap_unordered(atom_count_p, dataset),
            total=len(dataset),
            miniters=1000
        )
            if x is not None
        ]
    dataset = []
    dataset_cano = pd.DataFrame(dataset_cano, columns=['counts'])
    dataset_cano.to_csv(file_set+'_atom_counts.csv', index=None)
    # dataset_cano = dataset_cano.drop_duplicates('counts')
    # dataset_cano.to_csv(file_set+'.cano', index=None)


def main(config):
    #model = config.model
    #epoch = config.epoch
    sample_id = config.sample_id
    n_jobs = config.n_jobs
    path = config.path
    in_file = config.input
    output = config.output
    global dataset

    # file_set = file_collection(model, epoch)
    for sample_id in range(0,1,1):
        file_set = in_file  #'reinvent_1b_cano_diff.csv' # %(sample_id)
        dataset = open_file(file_set, path)
        # dataset = ['C12C3C4C5C4C4C(C4C13)C25',
        #            'C1C2C1C1C3C4CC4C2C13',
        #            'C1C2C1C1C3C4CC4C1C23',
        #            'C1C2C1C1C3C4CC(C24)C13',
        #            'C1C2C1C1C3CC4C2C4C31',
        #            'C1C2']
        counts(n_jobs, output, file_set)


if __name__ == '__main__':
    parser = get_parser()
    config, unknown = parser.parse_known_args()
    dataset = []
    main(config)
