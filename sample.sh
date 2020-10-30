# Tips 
model=$1
train_id=$2
epoch=$3
sample_batch=$4
gpu=$5
cd $model'_train_'$2
python -u  ../../moses_Apr08/scripts/sample.py $model \
       --model_load ./$1_model_s_$3.pt \
       --config_load ./$1_config \
       --vocab_load ./$1_voc \
       --n_samples  250000000 \
       --gen_save  $1_model_s_$3_$4.csv \
       --n_batch 20000  \
       --seed 100  \
       --device cuda:$5


