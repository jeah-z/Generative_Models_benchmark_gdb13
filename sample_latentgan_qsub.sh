for i in {10401..10500}
do
sh sample_latentgan.sh 0 $i  >sample_latentgan_0_1990_${i}_a.out  2>&1 
done
