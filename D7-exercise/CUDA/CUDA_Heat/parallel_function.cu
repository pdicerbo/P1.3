typedef float MYFLOAT;

#define PI 3.14159265359

__device__ MYFLOAT gpu_ix2x(int ix, int nx, MYFLOAT lx){
  return ((ix) - nx / 2.0)*lx/nx;
}

__device__ MYFLOAT gpu_iy2y(int iy, int ny, MYFLOAT ly){ 
  return ((iy) - ny / 2.0)*ly/ny;
}


/* function that initialize the temperature matrix on GPU. 
the initialization is done taking into account the presence of the ghosts cell
at the neighbour of the matrix. In this way, if SIZE is the size of the matrix without 
taking into accout the presence of the ghosts cell, to do the correct initialization
you must call this kernel with a number of threads and a number of blocks such that
NBLOCKS * THREADS_PER_BLOCK = SIZE */

__global__ void parallel_init(MYFLOAT *temp, int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay){ 

  int ix, iy, offset;
  MYFLOAT x, y;
  
  ix = threadIdx.x + blockIdx.x * blockDim.x;
  iy = threadIdx.y + blockIdx.y * blockDim.y;
  offset = ix + (iy + 1) * ( blockDim.x * gridDim.x + 2) + 1;
  
  x = gpu_ix2x(ix, nx, lx);
  y = gpu_iy2y(iy, ny, ly);
  
  temp[offset] = 1.0 / (2. * PI * sigmax * sigmay) * exp(-(x - x0)*(x - x0) / (2.0*(sigmax * sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay * sigmay)) );
}

__global__ void update_up_down(int nx, int ny, MYFLOAT *temp){

  int ix;
  
  ix = threadIdx.x + blockIdx.x * blockDim.x + 1;

  temp[ix] = temp[ix + blockDim.x * gridDim.x + 2];
  temp[(ny + 1) * (nx + 2) + ix] = temp[ny * (nx + 2) + ix];
}

__global__ void update_left_right(int nx, int ny, MYFLOAT *temp){

  int ix;
  
  ix = threadIdx.x + blockIdx.x * blockDim.x + 1;

  temp[ix * (nx + 2)] = temp[ix * (nx + 2) + 1];
  temp[ix * (nx + 2) + nx + 1] = temp[ix * (nx + 2) + nx]; //temp[ix * (nx + 2) + nx + 1];
}

