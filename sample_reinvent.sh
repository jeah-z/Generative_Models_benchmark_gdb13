#! /usr/bin/bash
# conda activate mosei
# mkdir reinvent_train_0_0
# cd reinvent_train_0_0
reinvent_path=/mnt/home/zhangjie/Mose_related/reinvent-randomized
nohup python $reinvent_path/sample_from_model.py -m ./models/model.trained.22 -n 1000000000 -o reinvent_1b_30_epoch.smi -b 12800 > sample_1b.out 2>&1 &
