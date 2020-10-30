#! /usr/bin/bash
# conda activate mosei
mkdir reinvent_train_Sep
cd reinvent_train_Sep
reinvent_path=/mnt/home/zhangjie/Mose_related/reinvent-randomized
python $reinvent_path/create_model.py -i /mnt/home/zhangjie/Mose_related/dataset_gdb13/$1/train.smi -o ./gdb_$1.ckpt  > create_$1.screen 2>&1 
nohup python $reinvent_path/train_model.py -i ./gdb_$1.ckpt -o ./Prior_$1/model -s /mnt/home/zhangjie/Mose_related/dataset_gdb13/$1/train.smi -e 1000 --lrm ada --csn 0  > train_$1.screen 2>&1 &
