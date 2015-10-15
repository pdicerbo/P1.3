#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 10

void print_lines(double*, int, int);

int main(int argc, char** argv){

  int MyRank, NPE, i, j, trasl, tmp_idx;
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

  if(rest != 0 && MyRank < rest){
    block += 1;
    trasl = block * MyRank;
  }
  else if(rest != 0)
    trasl = (block + 1) * rest + block * (MyRank - rest);
  else
    trasl = block * MyRank;

  tmp_buf = (double*)malloc(SIZE * block * sizeof(double));

  for(i = 0; i < block; i++){
    tmp_idx = i*SIZE;
    for(j = 0; j < SIZE; j++){
      if(j == trasl){
	tmp_buf[j + tmp_idx] = 1;
      }
      else
	tmp_buf[j + tmp_idx] = 0;
    }
    trasl += 1;
  }

  if(MyRank == 0){

#ifdef DEBUG

    if(SIZE <= 10){
      /* print matrix */
      printf("\n=============================\n");
      printf("\tPRINTING MATRIX:\n\n");
      
      print_lines(tmp_buf, SIZE, block);
      
      for(process = 1; process < NPE; process++){
	if( (rest != 0 && process < rest) || rest == 0 ){
	  MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	  print_lines(tmp_buf, SIZE, block);
	}
	else{
	  MPI_Recv(tmp_buf, SIZE*(block - 1), MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	  print_lines(tmp_buf, SIZE, block - 1);
	}
      }
      printf("\n=============================\n");
    }
    else{ 
      /* write in binary file */
      FILE *fp;
      fp = fopen("matrix.dat", "w");
      fwrite(tmp_buf, sizeof(double), SIZE*block, fp);
      for(process = 1; process < NPE; process++){
	if((rest != 0 && process < rest) || rest == 0){
	  MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	  fwrite(tmp_buf, sizeof(double), SIZE*block, fp);
	}
	else{
	  MPI_Recv(tmp_buf, SIZE*(block - 1), MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	  fwrite(tmp_buf, sizeof(double), SIZE*(block-1), fp);
	}
      }
    }

#else

    double* def_array;
    int *rec_count, *displs; 
    def_array = (double*)malloc(SIZE*SIZE*sizeof(double));
    rec_count = (int*)malloc(NPE*sizeof(int));
    displs = (int*)malloc(NPE*sizeof(int));

    rec_count[0] = SIZE * block;
    displs[0] = 0;

    for(process = 1; process < NPE; process++){
      if((rest != 0 && process < rest) || rest == 0){
	rec_count[process] = SIZE * block;
	displs[process] = SIZE * block + displs[process - 1];
      }
      else{
	rec_count[process] = SIZE * (block - 1);
	displs[process] = SIZE * (block - 1) + displs[process - 1];
      }
    }

    MPI_Gatherv(tmp_buf, SIZE*block, MPI_DOUBLE, def_array, rec_count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("\n\tDone\n");
    /* print_lines(def_array, SIZE, SIZE); */

#endif /* DEBUG */

  }

  else{

#ifdef DEBUG

    MPI_Send(tmp_buf, SIZE*block, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

#else

    int *rec_count, *displs;
    double *def_array;

    MPI_Gatherv(tmp_buf, SIZE*block, MPI_DOUBLE, def_array, rec_count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif /* DEBUG */
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
