#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 10

int main(int argc, char** argv){

  int MyRank, NPE, left, right, sum, data_rec, data_send, j;
  int tag = 42;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NPE);
  
  right = (MyRank + 1) % NPE;
  left = (MyRank + NPE - 1) % NPE;

  sum = MyRank;
  data_send = MyRank;
  j = 0;

  while(j < NPE - 1){
    MPI_Sendrecv(&data_send, 1, MPI_INT, right, tag, &data_rec, 1, MPI_INT, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    sum += data_rec;
    data_send = data_rec;
    j++;
  }
  
  /* printf("\n\tMyRank is %d of %d; my sum is %d; my left is %d; my right is %d;\n", MyRank, NPE, sum, left, right); */

  double *arr_rec, *arr_sum;
  int k;

  MPI_Request first_request, request;
  MPI_Status status;

  arr_rec  = (double*)malloc(SIZE * SIZE * sizeof(double));
  arr_sum  = (double*)malloc(SIZE * SIZE * sizeof(double));

  for(j = 0; j < SIZE; j++)
    arr_sum[j]  = MyRank;

  MPI_Isend(arr_sum, SIZE, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &first_request);

  j = 0;

  while(j < NPE - 1){
    MPI_Recv(arr_rec, SIZE, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &status);
    MPI_Isend(arr_rec, SIZE, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &request);

    MPI_Wait(&first_request, &status);

    for(k = 0; k < SIZE; k++)
      arr_sum[k] += arr_rec[k];

    j++;
  }

  printf("\n\tMyRank is %d of %d; my arr_sum[%d] = %lg;\n", MyRank, NPE, MyRank, arr_sum[MyRank]);

  MPI_Finalize();

  return 0;
}
