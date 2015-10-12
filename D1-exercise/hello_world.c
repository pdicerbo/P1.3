#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){
  
  int size, iam;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("\tI am %d of %d. Hello World!\n", iam, size);
  
  MPI_Finalize();

  return 0;
}
