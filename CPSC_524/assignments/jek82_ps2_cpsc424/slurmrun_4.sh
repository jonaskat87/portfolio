#!/bin/bash

#SBATCH --partition=day
#SBATCH --reservation=cpsc424
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=7G
#SBATCH --job-name=MandelbrotFinalTest
#SBATCH --output=%x-%j.out

# Load Required Modules

module load intel
module list

# Set side length, h, and number of iterations, N. Save into tmp.txt.
h=0.001
N=25000
echo $N > tmp.txt
echo $h >> tmp.txt

# Task 4

make clean
echo "make ts-random-2"
make ts-random-2
# Same as mandomp-collapse but modified to work with drand-ts-2.o
echo "make mandomp-collapse-ts"
make mandomp-collapse-ts

echo ""
echo ""
echo "OpenMP version with second version of threadsafe drand"
echo ""
export OMP_NUM_THREADS=24
export OMP_SCHEDULE="guided"
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE
for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandomp-collapse-ts
    echo ""
done

export OMP_SCHEDULE="static,1"
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE
for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandomp-collapse-ts
    echo ""
done
make clean
