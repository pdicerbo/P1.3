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
  
  int size, iam, imesg, i, nproc;
  double t_start, t_end, ndata, band;
  
  MPI_Init(&argc, &argv);

  nproc = atoi(argv[1]);

  double* test = (double*)malloc(SIZE*sizeof(double));

  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(iam == 0){
    /* initialization is done only by process 0 */
    for(i = 0; i < SIZE; i++){
      test[i] = ((double) i) * 0.5;
    }
  }

  /* imesg = iam * 10; */
  /* printf("\tBefore BCast operation, I am %d of %d. imesg is %d\n", iam, size, imesg); */
  
  /* MPI_Bcast(&imesg, 1, MPI_INT, 8, MPI_COMM_WORLD); */

  /* printf("\tAfter BCast operation, I am %d of %d. imesg is %d\n", iam, size, imesg); */

  if(iam == 0){
    t_start = seconds();
  }

  MPI_Bcast(test, SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(iam == 0){
    t_end = seconds();

    printf("\n================================================\n");
    /* printf("\tBANDWITH:\n"); */
    printf("\tArray size: %d; bytes / double = 8;\n\tnumber of procs: %d\n", SIZE, nproc);

    /* #data = SIZE * (byte / double) * #procs / (bytes -> Gbytes) */
    ndata = ((double)(SIZE))/((double)(1024 * 1024 *1024));
    /* this is needed to avoid int overflow */
    ndata *= 8 * nproc;

    printf("\tTotal data transfered: %lg GB", ndata);
    printf("\n\tTotal time: %lg s", (t_end - t_start));
    band = ndata / (t_end - t_start);
    printf("\n\tMesured Bandwidth: %lg GB/s\n", band);
    printf("\n================================================\n");
  }

  /* if(iam == 4){ */
  /*   printf("\n\t I am process %d!\n\tPrint vector:\n", iam); */
  /*   for(i = 0; i < SIZE; i++){ */
  /*     printf("%lg\t", test[i]); */
  /*   } */
  /* } */

  MPI_Finalize();

  return 0;
}
