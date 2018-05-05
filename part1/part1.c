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

static float smooth(float *input_arr, float *output_arr, int offset, int num_elements, int n, float a, float b, float c, int world_rank) {
  long int i, j = 0;
  long int n2 = n+2;
  float sum1, sum2, sum3;
  int row_offset = (offset / n);
  int col_offset = (offset % n);
  int output_index = 0;
  for (i = row_offset; i <=n; i++) {
    for (j = 1; j <=n; j++) {
      if (i == row_offset && j == 1 && col_offset > 1) {
        j = col_offset;
      }
      sum1  =  *(input_arr + (i-1)*n2 + j-1) + *(input_arr + (i-1)*n2 + j+1) + *(input_arr + (i+1)*n2 + j-1) + *(input_arr + (i+1)*n2 + j+1);
      sum1 = sum1*a;
      sum2 = *(input_arr + (i-1)*n2 + j) + *(input_arr + (i+1)*n2 + j) + *(input_arr + (i)*n2 + j-1) + *(input_arr + (i)*n2 + j+1);
      sum2 = sum2*b;
      sum3 = *(input_arr + (i)*n2 + j);
      sum3 = sum3*c;
      if (output_index == num_elements) {
        i = n+1;
        j = n+1;
        break;
      }
      else {
        output_index = output_index+1;
        *(output_arr + output_index) = sum1+sum2+sum3;
      }
    }
  }
  return *output_arr;
}

static void count(float *arr, int num_elements, float threshold, int *counter) {
  long int i, j;
  long int local_counter = 0;
  for (j = 1; j <=num_elements; j++) {
    if (*(arr + j) < threshold && *(arr + j) > 0) { //Do not count 0 elements - uninitialized elements in y
      local_counter++;
    }
  }
  *counter = local_counter;
}

void main(int argc, char *argv[]) { 

  //Declare all our stack variables
  float  a  = 0.05, b = 0.1, c = 0.4, t = 0.1;
  double times[6], start, end;
  long int n = 16386;
  int countx = 0, county = 0;

  //Init MPI
  int world_size;
  int world_rank;
  int irec = 5;
  MPI_Init(&argc,&argv);
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank); 

  //Init our seed for random number generation
  srand(time(NULL));

  //create x,y
  start = omp_get_wtime();
  float *x  = (float *)malloc( (n+2) * (n+2) * sizeof(float));
  end   = omp_get_wtime();
  times[0]  = end-start;

  start = omp_get_wtime();
  float *y  = (float *)malloc( (n+2) * (n+2) * sizeof(float));
  end   = omp_get_wtime();
  times[1]  = end-start;
  
  if (world_rank == 0) {
    start = omp_get_wtime();
    initialize(x, n);
    end   = omp_get_wtime();
    times[2]  = end-start;
  }
  //Send arr 'x' to all procs
  MPI_Bcast(x, (n+2)*(n+2), MPI_FLOAT, 0, MPI_COMM_WORLD);

  long int elements_per_proc = (n)*(n)/(world_size); //distribute uniformly among procs 
  float *recv_x       = (float *)malloc( elements_per_proc * sizeof(float));
  float *recv_y       = (float *)malloc( elements_per_proc * sizeof(float));

  MPI_Barrier(MPI_COMM_WORLD);
  start = omp_get_wtime();
  *recv_y = smooth(x, recv_y, elements_per_proc*world_rank, elements_per_proc, n, a, b, c, world_rank);
  MPI_Barrier(MPI_COMM_WORLD);
  end   = omp_get_wtime();
  times[3]  = end-start;

  // Scatter the random numbers to all processes
  MPI_Scatter(x+n+2, elements_per_proc, MPI_FLOAT, recv_x,
               elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  start = omp_get_wtime();
  // Compute the average of your subset
  count(recv_x, elements_per_proc, t, &countx);
  MPI_Barrier(MPI_COMM_WORLD);
  end   = omp_get_wtime();
  times[5]  = end-start;

  // Compute the average of your subset
  MPI_Barrier(MPI_COMM_WORLD);
  start = omp_get_wtime();
  count(recv_y, elements_per_proc, t, &county);
  MPI_Barrier(MPI_COMM_WORLD);
  end   = omp_get_wtime();
  times[4]  = end-start;

  // Gather all partial averages down to the root process
  int totalx;
  MPI_Reduce(&countx, &totalx, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  int totaly;
  MPI_Reduce(&county, &totaly, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Finalize();

  //Print info
  if (world_rank == 0) {
    fprintf(stdout, "Execution summary.");
    fprintf(stdout, "\n");
    fprintf(stdout, "Num elements in an array col                   : %i\n", n);
    fprintf(stdout, "Num inner elements in an array col             : %i\n", (n-2));
    fprintf(stdout, "Num elements in an array                       : %i\n", n * n);
    fprintf(stdout, "Num inner elements in an array                 : %i\n", (n-2) * (n-2));
    fprintf(stdout, "Num bytes used by an array (GB)                : %f\n", (double) (n * n * sizeof(float)) / 1000000000 ) ;
    fprintf(stdout, "Num elements below threshold in array y        : %i\n", totaly);
    fprintf(stdout, "Fraction of elements below threshold in array y: %f\n", (double) totaly / (double) (n*n));
    fprintf(stdout, "Num elements below threshold in array x        : %i\n", totalx);
    fprintf(stdout, "Fraction of elements below threshold in array x: %f\n", (double) totalx / (double) (n*n));
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
}
