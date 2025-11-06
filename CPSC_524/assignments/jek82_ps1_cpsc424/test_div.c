/* Benchmark division, in this case by dividing pi by increasing indices.*/

# include <stdio.h>
# include <stdlib.h> // for clearing screen, user input, and abs
# include "timing.h" // using provided timing routine
# include <math.h> // to get pi

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int main() {
  int N; // number of divisions
  double wcTime_i, wcTime_f, cpuTime, div; // start wallclock time, end wallclock time, cpu time, and division value, respectively
    
  // clear screen just in case for input
  system("clear");
  // take user input for N
  scanf("%d", &N);
  
  timing(&wcTime_i, &cpuTime); // start times for routine
  for (int i = 1; i < N; ++i) {
    div = M_PI / i;
  }
  timing(&wcTime_f, &cpuTime); // end times for routine

  // print value of div (for fun)
  printf("div=%.15lf.\n", div);
  // print wallclock time elapsed
  double wcTime = wcTime_f - wcTime_i;
  printf("Wallclock time elapsed for N divisions was %lf s.\n", wcTime);
  // print MFlops
  /* (note on how MFlops are computed in this case)
  There are N total floating-point divisions that have to be done, giving us
  N total floating-point operations. Hence, we get the Flops by dividing the
  result by the wallclock time, or wcTime. Finally, to turn Flops into 
  MFlops, we divide the result by 10^6=1000000. */
  printf("MFlops needed for numerical integration was %lf MFlops.\n", (double)N / wcTime / 1000000.0);
  
  return 0;
}
