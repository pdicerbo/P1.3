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

  double *send_buf, *rec_buf;
  double t_start, t_end;

  MPI_Request request[2];
  MPI_Status status[2];

  MPI_Init(&argc, &argv);

  send_buf = (double*)malloc(SIZE*sizeof(double));
  rec_buf  = (double*)malloc(SIZE*sizeof(double));


  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if(size != 2){
    printf("\n\tThis program must be executed with -np 2\n");
    exit(0);
  }

  /* OPTIMIZED VERSION: */
  /* INITIALIZATION */
  for(j = 0; j < SIZE; j++)
    send_buf[j] = ((double) j - SIZE)/SIZE;

  /* EXCHANGE DATA */
  t_start = seconds();

  MPI_Isend(send_buf, SIZE, MPI_DOUBLE, 1 - iam, tag, MPI_COMM_WORLD, &request[0]);
  MPI_Irecv(rec_buf, SIZE, MPI_DOUBLE, 1 - iam, tag, MPI_COMM_WORLD, &request[1]);
  MPI_Waitall(2, request, status);
  t_end = seconds();
  printf("\n===============================================\n");
  printf("\tRESULTS FROM PROCESS %d;\n", iam);
  printf("\tsend_buf[0] = %lg;\trec_buf[0] = %lg\n", send_buf[0], rec_buf[0]);
  printf("\ttotal time: %lg", (t_end - t_start));
  printf("\n===============================================\n");

  MPI_Finalize();

  return 0;
}
