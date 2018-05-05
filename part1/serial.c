#include <stdio.h>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

static void initialize(float *empty_arr, int n) {
  float random_float = 0.00;
  long int i, j    = 0;
  long int n2 = n+2;
  for (i = 0; i < n2; i++)
    for (j = 0; j < n2; j++) {
      random_float =(float) random()/RAND_MAX;
      *(empty_arr + i*n2 + j) = random_float;
    }
}
static float smooth(float *input_arr, float *output_arr, int n, float a, float b, float c) {
  long int i, j = 0;
  long int n2 = n+2;
  float sum1, sum2, sum3;
  for (i = 1; i <=n; i++) {
    for (j = 1; j <=n; j++) {
      sum1  =  *(input_arr + (i-1)*n2 + j-1) + *(input_arr + (i-1)*n2 + j+1) + *(input_arr + (i+1)*n2 + j-1) + *(input_arr + (i+1)*n2 + j+1);
      sum1 = sum1*a;
      sum2 = *(input_arr + (i-1)*n2 + j) + *(input_arr + (i+1)*n2 + j) + *(input_arr + (i)*n2 + j-1) + *(input_arr + (i)*n2 + j+1);
      sum2 = sum2*b;
      sum3 = *(input_arr + (i)*n2 + j);
      sum3 = sum3*c;
      *(output_arr + i*n2 + j) = sum1+sum2+sum3;
    }
  }
  return *output_arr;
}

static void count(float *arr, int n, float threshold, int *counter) {
  long int i, j;
  int local_counter = 0;
  for (i = 1;   i <=n; i++) {
    for (j = 1; j <=n; j++) {
      if (*(arr + i*(n+2) + j) < threshold) {
        local_counter++;
      }
    }
  }
  *counter = local_counter;
}

void main() { 

  //Declare all our stack variables
  float  a  = 0.05, b = 0.1, c = 0.4, t = 0.1;
  double times[6], start, end;
  int county, countx;
  long int n  = 16386;

  //Init our seed for random number generation
  srand(time(NULL));

  start = omp_get_wtime();
  float *x  = (float *)malloc( (n+2) * (n+2) * sizeof(float));
  end   = omp_get_wtime();
  times[0]  = end-start;

  start = omp_get_wtime();
  float *y  = (float *)malloc( (n+2) * (n+2) * sizeof(float));
  end   = omp_get_wtime();
  times[1]  = end-start;

  start = omp_get_wtime();
  initialize(x, n);
  end   = omp_get_wtime();
  times[2]  = end-start;

  start = omp_get_wtime();
  *y = smooth(x, y, n, a, b, c);
  end   = omp_get_wtime();
  times[3]  = end-start;

  start = omp_get_wtime();
  count(y, n, t, &county);
  end   = omp_get_wtime();
  times[4]  = end-start;

  start = omp_get_wtime();
  count(x, n, t, &countx);
  end   = omp_get_wtime();
  times[5]  = end-start;

  //Print info
  fprintf(stdout, "Execution summary.");
  fprintf(stdout, "\n");
  fprintf(stdout, "Num elements in an array col                   : %i\n", n);
  fprintf(stdout, "Num inner elements in an array col             : %i\n", (n-2));
  fprintf(stdout, "Num elements in an array                       : %i\n", n * n);
  fprintf(stdout, "Num inner elements in an array                 : %i\n", (n-2) * (n-2));
  fprintf(stdout, "Num bytes used by an array (GB)                : %f\n", (double) (n * n * sizeof(float)) / 1000000000 ) ;
  fprintf(stdout, "Num elements below threshold in array y        : %i\n", county);
  fprintf(stdout, "Fraction of elements below threshold in array y: %f\n", (double) county / (double) (n*n));
  fprintf(stdout, "Num elements below threshold in array x        : %i\n", countx);
  fprintf(stdout, "Fraction of elements below threshold in array x: %f\n", (double) countx / (double) (n*n));
  fprintf(stdout, "Threshold value:                               : %f\n",   t);

  //Timing each action
  fprintf(stdout, "\n");
  fprintf(stdout, "Timing summary (Seconds).");
  fprintf(stdout, "\n");
  fprintf(stdout, "WALL: Alloc-X :  %.6f\n", times[0]);
  fprintf(stdout, "WALL: Alloc-Y :  %.6f\n", times[1]);
  fprintf(stdout, "WALL: Init-X  :  %.6f\n", times[2]);
  fprintf(stdout, "WALL: Smooth-Y:  %.6f\n", times[3]);
  fprintf(stdout, "WALL: Count-Y :  %.6f\n", times[4]);
  fprintf(stdout, "WALL: Count-X :  %.6f\n", times[5]);

  //OMP Information.
  fprintf(stdout, "\n");
  fprintf(stdout, "OMP information.");
  fprintf(stdout, "\n");
  fprintf(stdout,"Num threads: %i\n", omp_get_max_threads());

}
