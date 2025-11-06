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

# Task 3 - Part 1

make clean
# mandomp-tasks replaces loop parallelism with task parallelism using a dedicated task generator
echo "make ts-random"
make ts-random
echo "make mandomp-tasks"
make mandomp-tasks

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with tasks"
echo ""
for no_threads in 1 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic,250"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-tasks
	echo ""
    done
done

# Task 3 - Part 2

make clean
# mandomp-tasks-row uses 1 entire row per task with a dedicated task generator
echo "make ts-random"
make ts-random
echo "make mandomp-tasks-row"
make mandomp-tasks-row

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with row tasks"
echo ""
for no_threads in 1 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic,250"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-tasks-row
	echo ""
    done
done

# Task 3 - Part 3

make clean
echo "make ts-random"
make ts-random
# Same as mandomp-tasks but shares the task generation among threads
echo "make mandomp-tasks-shared"
make mandomp-tasks-shared

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with shared task creation"
echo ""
for no_threads in 1 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic,250"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-tasks-shared
	echo ""
    done
done

# Task 3 - Part 4

make clean
echo "make ts-random"
make ts-random
# Same as mandomp-tasks-row but shares the task generation among threads
echo "make mandomp-tasks-row-shared"
make mandomp-tasks-row-shared

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with row tasks and shared task creation"
echo ""
for no_threads in 1 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic,250"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-tasks-row-shared
	echo ""
    done
done
make clean
