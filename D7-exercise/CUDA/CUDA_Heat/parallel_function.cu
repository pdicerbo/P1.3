typedef float MYFLOAT;

#define PI 3.14159265359

__device__ MYFLOAT gpu_ix2x(int ix, int nx, MYFLOAT lx){
  return ((ix) - nx / 2.0)*lx/nx;
}

__device__ MYFLOAT gpu_iy2y(int iy, int ny, MYFLOAT ly){ 
  return ((iy) - ny / 2.0)*ly/ny;
}

__global__ void update_up_down(int nx, int ny, MYFLOAT *temp){

  int ix;
  
  ix = threadIdx.x + blockIdx.x * blockDim.x + 1;

  temp[ix] = temp[ix + blockDim.x * gridDim.x + 2];
  temp[(ny + 1) * (nx + 2) + ix] = temp[ny * (nx + 2) + ix];

  // ix = threadIdx.x + blockIdx.x * blockDim.x;

  // temp[ix] = temp[ix + blockDim.x * gridDim.x];
  // temp[(ny + 1) * (nx + 2) + ix] = temp[ny * (nx + 2) + ix];
}

__global__ void update_left_right(int nx, int ny, MYFLOAT *temp){

  int ix;
  
  ix = threadIdx.x + blockIdx.x * blockDim.x + 1;

  temp[ix * (nx + 2)] = temp[ix * (nx + 2) + 1];
  temp[ix * (nx + 2) + nx + 1] = temp[ix * (nx + 2) + nx];

  // ix = threadIdx.x + blockIdx.x * blockDim.x;

  // temp[ix * (nx + 2)] = temp[ix * (nx + 2) + 1];
  // temp[ix * (nx + 2) + nx + 1] = temp[ix * (nx + 2) + nx];
}

__global__ void efficient_parallel_init(MYFLOAT *temp, int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay){ 

  int ix = threadIdx.x + blockIdx.x * ( blockDim.x - 2 );
  int iy = threadIdx.y + blockIdx.y * ( blockDim.y - 2 );
  int offset = ix + iy * (nx + 2); 
  
  MYFLOAT x = gpu_ix2x(ix - 1, nx, lx);
  MYFLOAT y = gpu_iy2y(iy - 1, ny, ly);
  
  temp[offset] = 1.0 / (2. * PI * sigmax * sigmay) * exp(-(x - x0)*(x - x0) / (2.0*(sigmax * sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay * sigmay)) );
}


__global__ void efficient_parallel_evolve(int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT *temp, MYFLOAT *temp_new, MYFLOAT alpha){
    
  MYFLOAT dx, dy;
  int ix = threadIdx.x;
  int iy = threadIdx.y;

  int g_ix = threadIdx.x + blockIdx.x * ( blockDim.x - 2);
  int g_iy = threadIdx.y + blockIdx.y * ( blockDim.y - 2);
  int offset = g_ix + g_iy * (nx + 2);
  
  dx = lx/nx;
  dy = ly/ny;
  
  extern __shared__ MYFLOAT temp_shrd[];
  extern __shared__ MYFLOAT temp_new_shrd[];
  
  temp_shrd[ix + iy * blockDim.x] = temp[offset];
    
  __syncthreads();
  
  if(ix > 0 && ix < blockDim.x - 1 && iy > 0 && iy < blockDim.y - 1){
    temp_new_shrd[ix + iy * blockDim.x] = temp_shrd[ix + iy * blockDim.x] + 
      alpha * dt *( (temp_shrd[ix + ( iy + 1 ) * blockDim.x] + 
		     temp_shrd[ix + ( iy - 1 ) * blockDim.x] - 
		     2.0 * temp_shrd[ix + iy * blockDim.x] ) / (dy * dy) +
		    (temp_shrd[ix + 1 + iy * blockDim.x] + 
		     temp_shrd[ix - 1 + iy * blockDim.x] - 
		     2.0 * temp_shrd[ix + iy * blockDim.x] ) / (dx * dx) );


    temp_new[offset] = temp_new_shrd[ix + iy * blockDim.x];
  }
}