#include <stdio.h>
#include <stdlib.h>

#define N 10
#define THREADS 10
#define BLOCKS 1

__global__ void reversing(double* arr, double* res, int size){
  int index = threadIdx.x + blockIdx.x * blockDim.x;

  res[size - 1 - index] = arr[index];
}

int main(int argc, char** argv){

  int i;
  double* arr = (double*)malloc(N * sizeof(double));
  double *dev_arr, *res;

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
  

  return 0;
}
