#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 10

void print_lines(double*, int, int);
void std_out_print(double*, int, int, int, int);

int main(int argc, char** argv){

  int MyID, NPE, i, j, i_global, tmp_idx;
  int block, left, right;

  int n_iter = 0;
  int i_in, j_in, k_in, offset, A_offset; /* inner index */

  double *A, *B, *C;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
  MPI_Comm_size(MPI_COMM_WORLD, &NPE);

  right = (MyID + 1) % NPE;
  left = (MyID + NPE - 1) % NPE;

  if(NPE > SIZE){
    printf("\n\tThe number of processes is bigger then the size of the problem\n");
    printf("\tNPE = %d > %d = SIZE\n\texit\n", NPE, SIZE);
    exit(0);
  }
  else if(SIZE % NPE != 0){
    printf("\n\tNot able to multiply matrix with SIZE % NPE != 0\n");
  }

  block = SIZE / NPE;

  A = (double*)malloc(SIZE * block * sizeof(double));
  B = (double*)malloc(SIZE * block * sizeof(double));
  C = (double*)malloc(SIZE * block * sizeof(double));

  for(i = 0; i < block; i++){

    tmp_idx = i*SIZE;
    i_global = block * MyID + i;

    for(j = 0; j < SIZE; j++){
      A[j + tmp_idx] = i_global - j;
      C[j + tmp_idx] = 0.;

      /* B is the identity matrix */
      if(j == i_global){
	A[j + tmp_idx] = 1.;
	B[j + tmp_idx] = 1.;
      }
      else{
	A[j + tmp_idx] = 0.;
	B[j + tmp_idx] = 0.;
      }
    }
  }


  /* loop over block */
  for(j = 0; j < SIZE; j += block){

    offset = j; // * block;

    /* single block loop */
    for(i_in = 0; i_in < block; i_in++){

      tmp_idx = i_in * SIZE;
      A_offset = MyID * block;

      for(j_in = offset; j_in < offset + block; j_in++){
	for(k_in = 0; k_in < block; k_in++)
	  C[j_in + tmp_idx] += A[k_in + tmp_idx + A_offset] * B[j_in + k_in*SIZE];
      }
    }
  }

  std_out_print(A, MyID, block, 0., NPE); /* The 0. stay for "rest" in previous std_out_print version*/

  std_out_print(B, MyID, block, 0., NPE); /* The 0. stay for "rest" in previous std_out_print version*/

  std_out_print(C, MyID, block, 0., NPE); /* The 0. stay for "rest" in previous std_out_print version*/

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
