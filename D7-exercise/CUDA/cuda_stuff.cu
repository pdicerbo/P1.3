#include <stdio.h>
#include <stdlib.h>

#define N 10
#define THREADS 10
#define BLOCKS 1

__global__ void reversing(double* arr, double* res, int size){
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  res[size - 1 - index] = arr[index];
}

__global__ void transpose(double* mat, double* res, int ncol){
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  
  // global index associated to matrix
  int r = index / ncol;
  int c = index % ncol;

  // transposition
  res[c * ncol + r] = mat[index];
}

__global__ void mm_mult(double* A, double* B, double* C, int nr_A, int nc_A, int nc_B){

  int index = threadIdx.x + blockIdx.x * blockDim.x;
  
  // global index associated to matrix
  int r = index / nc_B;
  int c = index % nc_B;
  int k;
  double tmp = 0.;

  for(k = 0; k < nc_A; k++)
    tmp += A[k + r * nc_A] * B[c + k * nc_A];

  C[index] = tmp;
}

int main(int argc, char** argv){

  int i, j, i_tmp;
  double* arr = (double*)malloc(N * sizeof(double));
  double *dev_arr, *res;

  double* mat = (double*)malloc(N * N * sizeof(double));
  double *dev_mat, *res_mat;
  double *identity = (double*)malloc(N * N * sizeof(double));
  double *dev_id, *dev_prod;

  if(THREADS * BLOCKS != N){
    printf("\n\tTHREADS * BLOCK must be equal to N\n\tExit\n");
    return 0;
  }
  
  cudaMalloc(&dev_arr, N * sizeof(double));
  cudaMalloc(&res, N * sizeof(double));

  for(i = 0; i < N; i++)
    arr[i] = (double) i;

  cudaMemcpy(dev_arr, arr, N * sizeof(double), cudaMemcpyHostToDevice);

  reversing<<<THREADS, BLOCKS>>>(dev_arr, res, N);

  cudaMemcpy(arr, res, N * sizeof(double), cudaMemcpyDeviceToHost);

  for(i = 0; i < N; i++)
    printf("\t%lg", arr[i]);
  printf("\n");

  cudaFree(dev_arr);
  cudaFree(res);
  free(arr);

  printf("\n\tTRANSPOSE\n\n");

  for(i = 0; i < N; i++){
    i_tmp = i * N;
    for(j = 0; j < N; j++){
      mat[j + i_tmp] = j + i_tmp;
      if(i != j)
	identity[i_tmp + j] = 0.;
      else
	identity[i_tmp + j] = 1.;
    }
  }

  for(i = 0; i < N; i++){
    i_tmp = i * N;

    for(j = 0; j < N; j++)
      printf("\t%lg", mat[i_tmp + j]);

    printf("\n");
}

  cudaMalloc(&dev_mat, N * N * sizeof(double));
  cudaMalloc(&res_mat, N * N * sizeof(double));
  cudaMalloc(&dev_id,  N * N * sizeof(double));
  cudaMalloc(&dev_prod,N * N * sizeof(double));

  cudaMemcpy(dev_mat, mat, N * N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_id, identity, N * N * sizeof(double), cudaMemcpyHostToDevice);

  transpose<<<THREADS, THREADS>>>(dev_mat, res_mat, N);
  mm_mult<<<THREADS, THREADS>>>(res_mat, dev_id, dev_prod, N, N, N);

  cudaMemcpy(mat, dev_prod, N * N * sizeof(double), cudaMemcpyDeviceToHost);
  printf("\n\n");
  for(i = 0; i < N; i++){

    i_tmp = i * N;

    for(j = 0; j < N; j++)
      printf("\t%lg", mat[i_tmp + j]);
    printf("\n");
  }
  printf("\n");

  return 0;
}
