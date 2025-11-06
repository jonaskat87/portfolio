/* testing out the range of drand() */

# include <stdlib.h>
# include <stdio.h>
# include "drand.h" // using provided random number generator
# include <omp.h>

int main() {
  int N; // number of threads

  // clear screen for input
  system("clear");
  // take user input for N
  scanf("%d", &N);

  dsrand(12345); // set random seed

  #pragma omp parallel num_threads(N)
  {
    #pragma omp for
    for (int i = 0; i < 20000; ++i) {
      printf("%lf\n", drand());
    }
  }
  return 0;
}
