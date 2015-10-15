#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 9

void print_lines(double*, int, int);

int main(int argc, char** argv){

  int MyRank, NPE, i, j, trasl, loop_ctrl;
  int block, rest, process;
  int tag = 42;

  double *tmp_buf;

  MPI_Request request;
  MPI_Status status;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Comm_size(MPI_COMM_WORLD, &NPE);

  if(NPE > SIZE){
    printf("\n\tThe number of processes is bigger then the size of the problem\n");
    printf("\tNPE = %d > %d = SIZE\n\texit\n", NPE, SIZE);
    exit(0);
  }

  block = SIZE / NPE;
  rest = SIZE % NPE;

  MyRank == 0 ? (loop_ctrl = 0) : (loop_ctrl = 1);

  if(rest != 0 && MyRank < rest){
    block += 1;
    trasl = (block - loop_ctrl + 1) * MyRank;
  }
  else if(rest != 0){
    trasl = (block + 1) * rest + block * (MyRank - rest);// (block - loop_ctrl + 1) * MyRank;
    /* trasl = (block + rest - loop_ctrl) * MyRank; */
  }
  else
    trasl = (block - loop_ctrl + 1) * MyRank;


  printf("\n\tMyRank: %d; block = %d; rest = %d; trasl = %d;\n", MyRank, block, rest, trasl);

  tmp_buf = (double*)malloc(SIZE * block * sizeof(double));

  for(i = 0; i < block; i++){
    for(j = 0; j < SIZE; j++){
      if(j == trasl){
      /* if(j == (MyRank + trasl)){ */
	tmp_buf[j + i*SIZE] = 1;
      }
      else
	tmp_buf[j + i*SIZE] = 0;
    }
    trasl += 1;
  }

  if(MyRank == 0){
    printf("\n=============================\n");
    printf("\tPRINTING MATRIX:\n");

    print_lines(tmp_buf, SIZE, block);

    for(process = 1; process < NPE; process++){
      if(rest != 0 && process < rest){
	MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	print_lines(tmp_buf, SIZE, block);
      }
      else if(rest != 0 && process >= rest){
	MPI_Recv(tmp_buf, SIZE*(block - 1), MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	print_lines(tmp_buf, SIZE, block - 1);
      }
      else{
	MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	print_lines(tmp_buf, SIZE, block);
      }
    }
  }
  else{
    MPI_Send(tmp_buf, SIZE*block, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
 
  MPI_Finalize();

  return 0;
}

void print_lines(double* buffer, int cols, int rows){
  int i, j, n_i;
  for(i = 0; i < rows; i++){
    n_i = i * cols;
    for(j = 0; j < cols; j++){
      printf("\t%lg", buffer[j + n_i]);
    }
    printf("\n");
  }
}
