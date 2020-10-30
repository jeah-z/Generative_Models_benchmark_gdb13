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


def disrib(gdb):
    perct = []
    size = len(gdb)
    
    gdb_head = gdb[gdb['counts']>20000]
    head_size = len(gdb_head)
    perct.append(['>20000', head_size/size*100])

    gdb_head = gdb[(gdb['counts']>2000) & (gdb['counts']<=20000)]
    head_size = len(gdb_head)
    perct.append(['2000~20000', head_size/size*100])

    gdb_head = gdb[(gdb['counts']>200) & (gdb['counts']<=2000)]
    head_size = len(gdb_head)
    perct.append(['200~2000', head_size/size*100])

    gdb_head = gdb[(gdb['counts']>20) & (gdb['counts']<=200)]
    head_size = len(gdb_head)
    perct.append(['20~200', head_size/size*100])

    gdb_head = gdb[(gdb['counts']>2) & (gdb['counts']<=20)]
    head_size = len(gdb_head)
    perct.append(['2~20', head_size/size*100])

    gdb_head = gdb[(gdb['counts']<=2)]
    head_size = len(gdb_head)
    perct.append(['<2', head_size/size*100])
    

    return perct



    # gt03 = df.sort_values('gt03_error', ascending=False)
    # for i in range(0, 1000):
    #     head = int(i*size/1000)
    #     print("%s was finished."%(str(i/1000*100)))
    #     gdb_head = gdb.head(head)
    #     # print(gdb_head['SMILES'])
    #     # print(model['SMILES'])
    #     same_flag = gdb_head['SMILES'].isin(model['SMILES'])
    #     # print(same_flag)
    #     same_flag=np.array(same_flag)
    #     same_flag.astype(int)
    #     # print(same_flag)
    #     print(sum(same_flag)/(len(same_flag)+1))
    #     perct.append([i/1000., sum(same_flag)/(len(same_flag)+1)])
    # return perct

ifg=1
rs=1
# ifg_path = 'aae_1b_ifg.csv'
# rs_path= 'aae_1b_rs.csv'


if ifg >0:
    pd_gdb = pd.read_csv("/mnt/home/zhangjie/Mose_related/gdb13_ifg_atom_nonfiltered_group.csv", names=['SMILES','counts'])
    percentage = disrib(pd_gdb)
    op_pert = pd.DataFrame(percentage, columns=["range", "percentage"])
    op_pert.to_csv('ifg_distribution.csv')

if rs >0:
    pd_gdb = pd.read_csv("/mnt/home/zhangjie/Mose_related/gdb13_ringsystem_nonfiltered_group.csv", names=['SMILES','counts'])
    percentage = disrib(pd_gdb)
    op_pert = pd.DataFrame(percentage, columns=["range", "percentage"])
    op_pert.to_csv('RS_distribution.csv')


