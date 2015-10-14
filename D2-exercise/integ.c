#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define SIZE 250000000

double seconds(){
  /* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) */
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

int main(int argc, char** argv){
  
  int size, iam, i, j, n_int_per_sub;
  double sum, w, large_h, x, total_sum;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* initializations */
  sum = 0;

  /* size of sub interval in which each process perform the sum */
  large_h = 1. / size;  

  /* number of sub sub interval in which each interval is divided */
  n_int_per_sub = 1e5;  //1e8 / size;

  /* size of each sub sub interval */
  w = large_h / n_int_per_sub;

  for(j = 1; j < n_int_per_sub; j++){
    x = large_h * iam + w * (j - 0.5);
    sum += 4 / (1 + x * x);
  }

  MPI_Reduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(iam == 0){
    /* printf("\n\tPi value: %lg\n", total_sum); */
    total_sum = total_sum / ((double) n_int_per_sub * size);
    printf("\n\tPi value: %.16g\n", total_sum);
  }

  MPI_Finalize();

  return 0;
}
