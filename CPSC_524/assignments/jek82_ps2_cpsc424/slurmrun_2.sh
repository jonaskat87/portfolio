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

# Task 2 - Part 1a

echo ""
echo ""
echo "OpenMP version with original drand.c"
echo ""

make clean
# mandomp uses original drand
echo "make mandomp-critical"
make mandomp-critical
echo "make mandomp-atomic"
make mandomp-atomic
echo "make mandomp-reduction"
make mandomp-reduction

export OMP_NUM_THREADS=2
unset OMP_SCHEDULE
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE

for i in {1..3}
do
    echo "Run "$i":"
    echo "Using critical:"
    cat tmp.txt | ./mandomp-critical
    echo "Using atomic:"
    cat tmp.txt | ./mandomp-atomic
    echo "Using reduction:"
    cat tmp.txt | ./mandomp-reduction
    echo ""
done

# Task 2 - Part 1b

make clean
# mandomp-ts is same as mandomp-reduction but uses threadsafe drand
echo "make ts-random"
make ts-random
echo "make mandomp-ts"
make mandomp-ts

echo ""
echo ""
echo "OpenMP version with threadsafe drand"
echo ""
for no_threads in 1 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    unset OMP_SCHEDULE
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-ts
	echo ""
    done
done

# Task 2 - Part 2a

make clean
echo "make ts-random"
make ts-random
echo "make mandomp-loops"
make mandomp-loops

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with loops, Part 2a: using schedule(static, 1)"
echo ""
for no_threads in 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="static,1"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-loops
	echo ""
    done
done

# Task 2 - Part 2b

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with loops, Part 2b: using schedule(static,100)"
echo ""
for no_threads in 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="static,100"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-loops
	echo ""
    done
done

# Task 2 - Part 2c

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with loops, Part 2c: using schedule(dynamic)"
echo ""
for no_threads in 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-loops
	echo ""
    done
done

# Task 2 - Part 2d

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with loops, Part 2d: using schedule(dynamic,250)"
echo ""
for no_threads in 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="dynamic,250"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-loops
	echo ""
    done
done

# Task 2 - Part 2e

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with loops, Part 2e: using schedule(guided)"
echo ""
for no_threads in 2 4 12 24;
do
    export OMP_NUM_THREADS=$no_threads
    export OMP_SCHEDULE="guided"
    echo "Number of threads = " $OMP_NUM_THREADS
    echo "OMP_SCHEDULE = " $OMP_SCHEDULE
    for i in {1..3}
    do
	echo "Run "$i":"
	cat tmp.txt | ./mandomp-loops
	echo ""
    done
done

# Task 2 - Part 3

make clean
echo "make ts-random"
make ts-random
echo "make mandomp-collapse"
make mandomp-collapse

echo ""
echo ""
echo "Performance runs for thread-safe OpenMP with collapse"
echo ""
export OMP_NUM_THREADS=24
export OMP_SCHEDULE="static,100"
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE
for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandomp-collapse
    echo ""
done

export OMP_NUM_THREADS=24
export OMP_SCHEDULE="dynamic,250"
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE
for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandomp-collapse
    echo ""
done

export OMP_NUM_THREADS=24
export OMP_SCHEDULE="guided"
echo "Number of threads = " $OMP_NUM_THREADS
echo "OMP_SCHEDULE = " $OMP_SCHEDULE
for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandomp-collapse
    echo ""
done
make clean
