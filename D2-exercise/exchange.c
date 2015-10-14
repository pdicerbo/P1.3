#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define SIZE 1000000000

double seconds(){
  /* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) */
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

int main(int argc, char** argv){

  int iam, size, j;
  int tag = 42;

  double *first_arr, *second_arr;
  double t_start, t_end;

  MPI_Request send_request, rec_request;
  MPI_Status status;
  MPI_Init(&argc, &argv);

  first_arr = (double*)malloc(SIZE*sizeof(double));
  second_arr = (double*)malloc(SIZE*sizeof(double));


  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if(size != 2){
    printf("\n\tThis program must be executed with -np 2\n");
    exit(0);
  }

  /* OPTIMIZED VERSION: */
  /* INITIALIZATION */
  for(j = 0; j < SIZE; j++)
    first_arr[j] = ((double) j - SIZE)/SIZE;

  /* EXCHANGE DATA */
  t_start = seconds();

  MPI_Isend(first_arr, SIZE, MPI_DOUBLE, 1 - iam, tag, MPI_COMM_WORLD, &send_request);
  MPI_Recv(second_arr, SIZE, MPI_DOUBLE, 1 - iam, tag, MPI_COMM_WORLD, &status);

  t_end = seconds();
  printf("\n===============================================\n");
  printf("\tRESULTS FROM PROCESS %d;\n", iam);
  printf("\tfirst_arr[0] = %lg;\tsecond_arr[0] = %lg\n", first_arr[0], second_arr[0]);
  printf("\ttotal time: %lg", (t_end - t_start));
  printf("\n===============================================\n");

  MPI_Finalize();

  /* INITIALIZATION */
  /* if(iam == 0){ */
  /*   for(j = 0; j < SIZE; j++) */
  /*     first_arr[j] = ((double) j - SIZE)/SIZE; */
  /* } */
  /* else{ */
  /*   for(j = 0; j < SIZE; j++) */
  /*     second_arr[j] = ((double) j + SIZE)/SIZE; */
  /* } */

  /* /\* EXCHANGE DATA *\/ */
  /* if(iam == 0){ */
  /*   MPI_Isend(first_arr, SIZE, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &send_request); */
  /*   MPI_Recv(second_arr, SIZE, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &status); */
  /*   /\* MPI_Irecv(second_arr, SIZE, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &rec_request); *\/ */
  /* } */

  /* if(iam == 1){ */
  /*   t_start = seconds(); */
  /*   MPI_Isend(second_arr, SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &send_request); */
  /*   MPI_Recv(first_arr, SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status); */
  /*   /\* MPI_Irecv(first_arr, SIZE, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &rec_request); *\/ */
  /*   /\* MPI_Wait(&rec_request, &status); *\/ */
  /*   t_end = seconds(); */
  /*   printf("\n\tRESULTS FROM PROCESS %d;\n", iam); */
  /*   printf("\tfirst_arr[0] = %lg;\tsecond_arr[0] = %lg\n", first_arr[0], second_arr[0]); */
  /*   printf("\ttotal time: %lg\n=================================\n", (t_end - t_start)); */
  /* } */

  return 0;
}
