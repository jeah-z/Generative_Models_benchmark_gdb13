# Generative_Models_benchmark_gdb13

This repository includes the script used training, sampling and analyzing of generative models in the project "Comparative study of deep generative models on chemical space coverage"

![image](https://github.com/jeah-z/Generative_Models_benchmark_gdb13/blob/master/summary/images/Coverage%20of%20GDB-13%20from%201B%20Sampled%20Compounds.png)
Fig. 1 Coverage of GDB-13 from 1B Sampled Compounds

![image](https://github.com/jeah-z/Generative_Models_benchmark_gdb13/blob/master/summary/images/Distribution%20of%20ring%20systems%20and%20functional%20groups%20in%20GDB-13.png)
Fig. 2 Distribution%20of%20ring%20systems%20and%20functional%20groups%20in%20GDB-13

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

# Related github repository

- CharRNN, AAE, VAE, ORGAN models were retreived from Github repository [https://github.com/molecularsets/moses]
- REINVENT was retreived from Github repository [https://github.com/undeadpixel/reinvent-randomized]
- LatentGAN was retreived from Github repository [https://github.com/Dierme/latent-gan]
- GraphINVENT was retreived from Github repository [https://github.com/MolecularAI/GraphINVENT]