#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 10

void print_lines(double*, int, int);
void std_out_print(double*, int, int, int, int);

int main(int argc, char** argv){

  int MyRank, NPE, i, j, i_global, tmp_idx;
  int block, rest, offset;

  double *tmp_buf;

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
  offset = 0;

  if(rest != 0 && MyRank < rest)
    block += 1;
  else
    offset = rest;

  tmp_buf = (double*)malloc(SIZE * block * sizeof(double));

  for(i = 0; i < block; i++){

    tmp_idx = i*SIZE;
    i_global = block * MyRank + offset + i;

    for(j = 0; j < SIZE; j++){
      if(j == i_global){
	tmp_buf[j + tmp_idx] = 1;
      }
      else
	tmp_buf[j + tmp_idx] = 0;
    }
  }

  std_out_print(tmp_buf, MyRank, block, rest, NPE);

  MPI_Finalize();

  return 0;
}

void std_out_print(double* tmp_buf, int MyRank, int block, int rest, int NPE){

  int process;

  int tag = 42;
  MPI_Status status;

  if(MyRank == 0){

#ifdef DEBUG

    if(SIZE <= 10){
      /* print matrix */
      printf("\n=============================\n");
      printf("\tPRINTING MATRIX:\n\n");
      
      print_lines(tmp_buf, SIZE, block);
      
      for(process = 1; process < NPE; process++){
	if(rest != 0 && process == rest)
	  block--;
	MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	print_lines(tmp_buf, SIZE, block);
      }
      printf("\n=============================\n");
    }
    else{ /* SIZE > 10 */
      /* write in binary file */
      FILE *fp;
      fp = fopen("matrix.dat", "w");
      fwrite(tmp_buf, sizeof(double), SIZE*block, fp);

      for(process = 1; process < NPE; process++){
	if(rest != 0 && process == rest)
	  block--;
	MPI_Recv(tmp_buf, SIZE*block, MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
	fwrite(tmp_buf, sizeof(double), SIZE*block, fp);
      }
      printf("\n=============================\n");
      printf("\n\tData stored in matrix.dat\n");
      printf("\n=============================\n");
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
      if(rest != 0 && process == rest)
	block--;
      rec_count[process] = SIZE * block;
      displs[process] = rec_count[process - 1] + displs[process - 1];
    }

    block++;

    MPI_Gatherv(tmp_buf, SIZE*block, MPI_DOUBLE, def_array, rec_count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("\n=============================\n");
    printf("\n\tDone\n");
    /* print_lines(def_array, SIZE, SIZE); */
    printf("\n=============================\n");

    free(rec_count);
    free(displs);
    free(def_array);

#endif /* DEBUG */

  }

  else{ /* MyRank > 0 */

#ifdef DEBUG

    MPI_Send(tmp_buf, SIZE*block, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

#else

    int *rec_count, *displs;
    double *def_array;

    MPI_Gatherv(tmp_buf, SIZE*block, MPI_DOUBLE, def_array, rec_count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#endif /* DEBUG */
  }
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
