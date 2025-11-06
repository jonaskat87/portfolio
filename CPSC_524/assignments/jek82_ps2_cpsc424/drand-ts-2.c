#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <omp.h>

static uint64_t seed;
// set seed to threadprivate so each thread has its own copy
#pragma omp threadprivate(seed)

void dsrand(unsigned s)
{
  // no need for plus one here
  seed = s;
}

/* The static variable seed has already been (or should be) declared 
   privately within each separate thread using thread private. This means 
   that drand is called on the private copy of seed in each thread. */

double drand(void)
{
  #pragma omp atomic
  seed *= 6364136223846793005ULL;
  #pragma omp atomic
  seed += 1;
  return((double)(seed>>33)/(double)RAND_MAX);
}
