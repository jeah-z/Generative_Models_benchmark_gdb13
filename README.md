# Generative_Models_benchmark_gdb13

## Requirements:
* rdkit
* collections
* numpy 
* pandas
* functools
* multiprocessing
* tqdm
* pathlib
* argparse


## Usage:

To calculate the function groups of the SMILES in the file gdb13_split_0.csv.
```
python cal_ifg_atom.py -n 40 -i gdb13_split_0
```

To calculate the ringsystems of the SMILES in the file gdb13_split_0.csv.
```
python cal_ringsystem.py -n 40 -i gdb13_split_0
```
- -n   number of threads will be used to process the file.
- -i   name of the input file without '.csv'


To count the atoms in the compounds of train.smi.
```
python atom_counts.py -n 2 -i Dataset/train.smi
```
After calculation, all the total number of atoms in the compounds of input file will be saved into another file with a file suffix '_atom_counts.csv'. 