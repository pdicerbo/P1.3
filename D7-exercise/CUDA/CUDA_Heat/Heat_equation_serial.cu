/*
 *     Using Time Forward Space Centered finite difference method this code solves the heat equation
 *                            u,t = alpha* (u,xx + u,yy) 
 *
 *    Boundaries are flat: copying the value of the neighbour
 *    Created by G.P. Brandino for MHPC
 *    Last Revision: October 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#ifdef LINK_OMP
#include <omp.h>
#endif

#define PI 3.14159265359

typedef float MYFLOAT;

#ifdef __CUDA
#define BLOCKS 2
#define THREADS 5
__global__ void parallel_init(MYFLOAT*, int, int, MYFLOAT, MYFLOAT, MYFLOAT, MYFLOAT, MYFLOAT, MYFLOAT);
#endif

/* 
 * conversions from discrete to real coordinates
 */
MYFLOAT ix2x(int ix, int nx, MYFLOAT lx){
    return ((ix-1) - nx / 2.0)*lx/nx;
}

MYFLOAT iy2y(int iy, int ny, MYFLOAT ly){ 
    return ((iy-1) - ny / 2.0)*ly/ny;
}


/* Integral of the temperature field */
MYFLOAT integral(int nx, int ny, MYFLOAT* val, MYFLOAT lx, MYFLOAT ly){
	MYFLOAT sum=0.0;
	int ix,iy;
#ifdef LINK_OMP
#pragma omp parallel for private(ix) collapse(2) reduction (+:sum)
#endif
    	for(iy=1;iy<=ny;++iy)
		for(ix=1;ix<=nx;++ix){
          		sum+=val[((nx+2)*iy)+ix];
 		}

 	return(sum*lx/nx*ly/ny);
}


/* 
 * initialize the system with a gaussian temperature distribution. The center of the gaussian (x0,y0) as well as the amplitudes sigmax, sigmay are input parameters 
 */
int init(MYFLOAT *temp, int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay){ 
    int ix,iy;
    MYFLOAT x,y;
#ifdef LINK_OMP
#pragma omp parallel for private(x, y, ix) collapse(2)
#endif
    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
	    x=ix2x(ix, nx, lx);
	    y=iy2y(iy, ny, ly);
	    temp[((nx+2)*iy)+ix] = 1.0/(2.*PI*sigmax*sigmay)*exp(-(x-x0)*(x-x0)/(2.0*(sigmax*sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay*sigmay)) );
	}

    printf("--  Max (for gnuplot) -- %f\n", 1.0/(2.*PI*sigmax*sigmay) );
    return 0;
}

/*
 * save the temperature distribution
 * the ascii format is suitable for splot gnuplot function
 */
int save_gnuplot(FILE* fp, MYFLOAT *temp, int nx, int ny, MYFLOAT lx, MYFLOAT ly) {
    
    int ix,iy;
    // the double newline will be interpreted by gnuplot as a new data block
    fprintf(fp, "\n\n");
    for(iy=1;iy<=ny;++iy){		
	for(ix=1;ix<=nx;++ix)
	    fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy,ny,ly), temp[((nx+2)*iy)+ix]);
    }

    return 0;
}

/*
 * Updating boundaries. Boundaries are flat: copying the value of the neighbour
 */

void update_boundaries_FLAT(int nx, int ny, MYFLOAT *temp){

    int ix, iy;
#ifdef LINK_OMP
#pragma omp parallel
#endif

#ifdef LINK_OMP
#pragma omp for
#endif
    for(iy=0;iy<=ny+1;++iy){
	temp[(nx+2)*iy] = temp[((nx+2)*iy)+1];
	temp[((nx+2)*iy)+(nx+1)] = temp[((nx+2)*iy)+nx];
    }

#ifdef LINK_OMP
#pragma omp for
#endif
    for(ix=0;ix<=nx+1;++ix){
	temp[ix] = temp[(nx+2)+ix];
	temp[((nx+2)*(ny+1))+ix] = temp[((nx+2)*ny)+ix];
    }
}

/*
 * Perform time evolution using FTCS FD method
 */

void evolve(int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT *temp, MYFLOAT *temp_new, MYFLOAT alpha){
    
    MYFLOAT dx, dy;
    int ix, iy;

    dx = lx/nx;
    dy = ly/ny;

#ifdef LINK_OMP
#pragma omp parallel for private(ix) collapse(2)
#endif
    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
	    temp_new[((nx+2)*iy)+ix] = temp[((nx+2)*iy)+ix] + alpha*dt*
              ( (temp[((nx+2)*(iy+1))+ix] + temp[((nx+2)*(iy-1))+ix] -2.0* temp[((nx+2)*iy)+ix])/(dy*dy) + 
                (temp[((nx+2)*iy)+(ix+1)] + temp[((nx+2)*iy)+(ix-1)] -2.0* temp[((nx+2)*iy)+ix])/(dx*dx) );
	}

#ifdef LINK_OMP
#pragma omp parallel for private(ix) collapse(2)
#endif
    for(iy=0;iy<=ny+1;++iy)
	for(ix=0;ix<=nx+1;++ix)
	    temp[((nx+2)*iy)+ix] = temp_new[((nx+2)*iy)+ix];
}

double seconds(){
  /* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) */
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

int main(int argc, char* argv[]){

  int i, j, nx, ny, n_steps, nRows, nCols, frame;
  MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
  MYFLOAT *temp, *temp_new;
  MYFLOAT *dev_temp, *dev_temp_new;
  MYFLOAT norm, norm_ini, bound;
  FILE *fp;
  
#ifdef LINK_OMP
  FILE *agg;
#endif
  
  double t_start, t_end, t_sum = 0.;
  
  if (argc !=3) 
    {
      printf(" FTCS finite difference solution of the heat equation \n\n");
      printf(" Usage : %s n_steps dt \n", argv[0]);
      return 1;
    }
  
  // number of points in the x directions
  nx = 10; //4000;
  nCols= nx + 2;
  
  // number of points in the y directions
  ny = 10; //4000;
  nRows= ny + 2;
  
  // size of the system in the x direction
  lx = 2.0;
  
  // size of the system in the y direction
  ly = 2.0;
  
  //center of the initial gaussian distribution    
  x0 = 0.0;
  y0 = 0.0;
  
  // amplitude of the initial gaussian distribution 
  sigmax = 0.5;
  sigmay = 0.2;
  
  // alpha coefficient (diffusivity) of the heat equation 
  alpha = 1.0;
  
  // number of time steps and amplitude thereof. Command line argumentes 
  n_steps = atoi(argv[1]);
  dt = atof(argv[2]);
  
  // checking stability of the method 
  bound = 1.0/(2.0*alpha)*(1.0/(nx*nx/(lx*lx)+ny*ny/(ly*ly)));
  if (dt>bound  )
    printf("Warning, unstable regime, reduce dt \n");
  

#ifdef __CUDA

  temp = (MYFLOAT *)malloc(nRows * nCols * sizeof(MYFLOAT));

  cudaMalloc((void**)&dev_temp, nRows * nCols * sizeof(MYFLOAT));
  cudaMalloc((void**)&dev_temp_new, nRows * nCols * sizeof(MYFLOAT));

  dim3 blocks(BLOCKS, BLOCKS);
  dim3 threads(THREADS, THREADS);
  // dim3 blocks, threads;
  // blocks.x = BLOCKS;
  // blocks.y = BLOCKS;
  // threads.x = THREADS;
  // threads.y = THREADS;

  parallel_init<<<blocks, threads>>>(dev_temp, nx, ny, lx, ly, x0, y0, sigmax, sigmay);
  cudaMemcpy(temp, dev_temp, nRows * nCols * sizeof(MYFLOAT), cudaMemcpyDeviceToHost);

  for(i = 1; i <= ny; i++){
    for(j = 1; j <= nx; j++)
      printf("\t%lg", temp[(nx + 2) * i + j]);
    printf("\n");
  }

  temp_new = (MYFLOAT *)malloc(nRows * nCols * sizeof(MYFLOAT));
  printf("\n--------------------\n\n");
  // fill temperaature array with initial condition, imposing flat boundary conditions
  init(temp_new, nx, ny, lx, ly, x0, y0, sigmax, sigmay);
  for(i = 1; i <= ny; i++){
    for(j = 1; j <= nx; j++)
      printf("\t%lg", temp_new[(nx + 2) * i + j]);
    printf("\n");
  }
  exit(0);

#endif

  temp = (MYFLOAT *)malloc(nRows * nCols * sizeof(MYFLOAT));
  temp_new = (MYFLOAT *)malloc(nRows * nCols * sizeof(MYFLOAT));

  // fill temperaature array with initial condition, imposing flat boundary conditions
  init(temp, nx, ny, lx, ly, x0, y0, sigmax, sigmay);
  update_boundaries_FLAT(nx, ny,temp);
  
  // calculating initial norm (~ total energy)
  norm_ini=integral(nx, ny, temp, lx, ly);
  printf(" Initial integral value: %f\n", norm_ini);
  
  fp = fopen("heat_diffusion.dat", "w");

  frame=n_steps;
  printf(" Starting time evolution... \n\n ");
  t_start = seconds();
  for(i=1; i<=n_steps; ++i) {
    // saving a snapshot every 100 time steps
    if ( (i-1)%frame==0){
      t_end = seconds();
      save_gnuplot(fp, temp, nx, ny, lx, ly);
      t_sum += t_end - t_start;
      t_start = seconds();
    }
    // performing TFSC-FD step and updating boundaries
    evolve(nx, ny, lx, ly, dt, temp, temp_new, alpha);
    update_boundaries_FLAT(nx, ny, temp);
  }
  t_end = seconds();
  t_sum += t_end - t_start;
  
  // checking final norm (~total energy)
  norm=integral(nx, ny, temp, lx, ly); 
  printf(" Integral end: %f\n",norm);
#ifdef LINK_OMP
  printf("\t-----------------------\n");
  printf("\tNumber of threads: %d\n", omp_get_max_threads());
#endif
  printf("\tTime used: %lg s\n\n", t_sum);
#ifdef LINK_OMP
  agg = fopen("aggiorn.dat", "a");
  fprintf(agg, "\nrun with %d threads. Time used = %lg s.", omp_get_max_threads(), t_sum);
  fclose(agg);
#endif
  fclose(fp); 
  
  return 0;
}

