#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#define SIZE 2400

void print_lines(double*, int, int);
void std_out_print(double*, int, int, int, int);
void partial_prod(double*, double*, double*, int, int);
double seconds();

int main(int argc, char** argv){

  int MyID, NPE, i, j, k, i_global, tmp_idx;
  int block, prev, next;

  int shift;
  int tag = 42;

  double *A, *B, *C, *B_rec, *B_tmp;
  double t_start, t_end, t_wait_start, t_wait_end, t_compute, total_wait = 0;

  MPI_Request first_request, request = MPI_REQUEST_NULL;
  MPI_Status status;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
  MPI_Comm_size(MPI_COMM_WORLD, &NPE);

  next = (MyID + 1) % NPE;
  prev = (MyID + NPE - 1) % NPE;

  shift = MyID;

  if(NPE > SIZE){
    printf("\n\tThe number of processes is bigger then the size of the problem\n");
    printf("\tNPE = %d > %d = SIZE\n\texit\n", NPE, SIZE);
    exit(0);
  }
  else if(SIZE % NPE != 0){
    printf("\n\tNot able to multiply matrix with SIZE % NPE != 0\n");
    exit(0);
  }

  block = SIZE / NPE;

  A = (double*)malloc(SIZE * block * sizeof(double));
  B = (double*)malloc(SIZE * block * sizeof(double));
  C = (double*)malloc(SIZE * block * sizeof(double));
  
  for(i = 0; i < block; i++){

    tmp_idx = i*SIZE;
    i_global = block * MyID + i;

    for(j = 0; j < SIZE; j++){

      /* A[j + tmp_idx] = i_global * SIZE + j + 1; */
      /* B[j + tmp_idx] = i_global * SIZE + j + 1; */

      /* timing initialization */
      A[j + tmp_idx] = (double) (i_global + SIZE - j);
      /* B[j + tmp_idx] = (double) (i_global - j); */
      C[j + tmp_idx] = 0.;

      /* B (for semplicity) is the identity matrix */
      if(j == i_global){
      	/* A[j + tmp_idx] = 1.; */
      	B[j + tmp_idx] = 1.;
      }
      else{
      	/* A[j + tmp_idx] = 0.; */
      	B[j + tmp_idx] = 0.;
      }
    }
  }

  partial_prod(A, B, C, block, shift);

  B_rec = (double*)malloc(SIZE * block * sizeof(double));
  B_tmp = (double*)malloc(SIZE * block * sizeof(double));


  if(MyID == 0)
    t_start = seconds();

  MPI_Isend(B, SIZE * block, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, &first_request);

  for(i = 1; i < NPE; i++){

    MPI_Recv(B_rec, SIZE * block, MPI_DOUBLE, next, tag, MPI_COMM_WORLD, &status);

    if(MyID == 0)
      t_wait_start = seconds();

    MPI_Wait(&request, &status);

    if(MyID == 0){
      t_wait_end = seconds();
      total_wait += t_wait_end - t_wait_start;
    }

    for(k = 0; k < block; k++){
      tmp_idx = k * SIZE;
      for(j = 0; j < SIZE; j++)
	B[j + tmp_idx] = B_rec[j + tmp_idx];
    }

    if(i < NPE - 1)
      MPI_Isend(B, SIZE * block, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, &request);
    
    shift = (shift + 1) % NPE;
    partial_prod(A, B_rec, C, block, shift);
  }

  if(MyID == 0){
    t_end = seconds();
    t_compute = t_end - t_start - total_wait;
    printf("\n\tMatrix size: %d;\tCommunication time: %lg;\tCompute time: %lg\n", SIZE, total_wait, t_compute);

    FILE *fp;
    fp = fopen("timing.dat", "a");
    
    fprintf(fp, "%d\t%d\t%lg\t%lg\n", NPE, SIZE, total_wait, t_compute);

    fclose(fp);

  }
  /* Printing result */
  /* The 0. stay for "rest" in previous std_out_print version*/

  /* std_out_print(A, MyID, block, 0., NPE); */
  /* std_out_print(B, MyID, block, 0., NPE); */
  /* std_out_print(C, MyID, block, 0., NPE); */

  MPI_Finalize();

  return 0;
}

void partial_prod(double* A, double* B, double* C, int block, int shift){

  int j, i_in, j_in, k_in, offset, A_offset, tmp_idx; /* inner-block index */

  /* loop over blocks */
  for(j = 0; j < SIZE; j += block){

    /* single block loop */
    for(i_in = 0; i_in < block; i_in++){

      tmp_idx = i_in * SIZE;
      A_offset = shift * block;

      for(j_in = j; j_in < j + block; j_in++)
	for(k_in = 0; k_in < block; k_in++)
	  C[j_in + tmp_idx] += A[k_in + tmp_idx + A_offset] * B[j_in + k_in*SIZE];
    }
  }  
}

void std_out_print(double* tmp_buf, int MyRank, int block, int rest, int NPE){

  int process;

  int tag = 42;
  MPI_Status status;

  if(MyRank == 0){

#ifdef DEBUG

    if(SIZE <= 40){
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

double seconds(){
  /* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) */
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}
