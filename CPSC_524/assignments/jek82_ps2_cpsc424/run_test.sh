# specify number of cores
cores=12
echo $cores | ./rand_test | awk '{if (NR>1) print}' | tr '\n' ' ' | tr ' ' ',' > rand_samples.txt
