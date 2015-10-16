#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char** argv){

  int MyRank, NPE, left, right, sum, data_rec, data_send, j;
  int tag = 42;

  int l_ctrl;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NPE);
  
  right = (MyRank + 1) % NPE;
  l_ctrl = 1 / (MyRank + 1); /* == 1 if MyRank==0; == 0 Otherwise */
  left = (MyRank - 1) * (1 - l_ctrl) + (NPE - 1) * l_ctrl;

  sum = MyRank;
  data_send = MyRank;
  j = 0;

  while(j < NPE - 1){
    MPI_Sendrecv(&data_send, 1, MPI_INT, right, tag, &data_rec, 1, MPI_INT, left, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    sum += data_rec;
    data_send = data_rec;
    j++;
  }
  
  printf("\n\tMyRank is %d of %d; my sum is %d; my left is %d; my right is %d;\n", MyRank, NPE, sum, left, right);

  MPI_Finalize();

  return 0;
}