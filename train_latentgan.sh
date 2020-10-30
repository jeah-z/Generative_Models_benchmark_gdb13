# Tips   train.sh aae/char_rnn/vae/organ    0/1/2 
model=latentgan
train_id=$1
mkdir $model'_train_'$1
cd $model'_train_'$1
python /mnt/home/zhangjie/Mose_related/latent-gan-master/run.py  \
       --smiles-file ../../dataset_gdb13/$1/train.csv \
       --storage-path ./ \
       --latent-file ./latentgan_$0 \
       --n-epochs 2000  \
       --sample-n 0

 

       

