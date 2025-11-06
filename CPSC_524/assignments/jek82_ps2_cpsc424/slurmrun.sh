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

# Task 1

make clean
echo "make mandseq"
make mandseq

echo ""
echo ""
echo "Serial version"

for i in {1..3}
do
    echo "Run "$i":"
    cat tmp.txt | ./mandseq
    echo ""
done

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
echo "make ts-random"
make ts-random
# mandomp-ts is same as mandomp-reduction but uses threadsafe drand
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
export OMP_SCHEDULE="static,10"
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

# Task 4

make clean
echo "make ts-random-2"
make ts-random-2
# Same as mandomp-collapse but modified to work with drand-ts-2.o
echo "make mandomp-collapse-ts"
make mandomp-collapse-ts

echo ""
echo ""
echo "OpenMP version with second version of threadsafe drand and collapse"
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
