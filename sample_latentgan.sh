# Tips   train.sh aae/char_rnn/vae/organ    0/1/2 
model=latentgan
train_id=$1
#mkdir $model'_train_'$1
cd $model'_train_'$1
python /mnt/home/zhangjie/Mose_related/latent-gan-master/sample.py  \
       --generator-path ./1990_generator.txt \
       -olf ./latentgan_model_s_1990_latent_$2_a \
       --decode-sampled True \
       --number-samples 500000  \
       -odsf ./latentgan_model_s_1990_$2_a.csv

 

       

