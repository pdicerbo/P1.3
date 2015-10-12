#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){
  
  int size, iam, imesg;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &iam);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  imesg = iam * 10;
  printf("\tBefore BCast operation, I am %d of %d. imesg is %d\n", iam, size, imesg);
  
  MPI_Bcast(&imesg, 1, MPI_INT, 8, MPI_COMM_WORLD);

  printf("\tAfter BCast operation, I am %d of %d. imesg is %d\n", iam, size, imesg);

  MPI_Finalize();

  return 0;
}
