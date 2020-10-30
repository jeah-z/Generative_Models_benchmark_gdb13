# Tips   train.sh aae/char_rnn/vae/organ    0/1/2 
model=$1
train_id=$2
mkdir $model'_train_'$2
cd $model'_train_'$2
python ../../moses_Apr08/scripts/train.py $model \
       --train_load ../../dataset_gdb13/$2/train.csv \
       --val_load ../../dataset_gdb13/$2/test.csv \
       --model_save ./$1_model_save \
       --save_frequency 1 \
       --config_save ./$1_config \
       --vocab_save ./$1_voc \
       --log_file ./$1_$2.log  \
       --device cuda:1
       
       

