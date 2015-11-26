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
#include <omp.h>

#define PI 3.14159265359

typedef double MYFLOAT;

#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>


double cclock()
/* Returns elepsed seconds past from the last call to timer rest */
{

  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}


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

    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
	    x=ix2x(ix, nx, lx);
	    y=iy2y(iy, ny, ly);
	    temp[((nx+2)*iy)+ix] = 1.0/(2.*PI*sigmax*sigmay)*exp(-(x-x0)*(x-x0)/(2.0*(sigmax*sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay*sigmay)) );
	}
    return 0;
}


/*
 * Updating boundaries. Boundaries are flat: copying the value of the neighbour
 */

void update_boundaries_FLAT(int nx, int ny, MYFLOAT *temp){

    int ix, iy;

    for(iy=0;iy<=ny+1;++iy){
	temp[(nx+2)*iy] = temp[((nx+2)*iy)+1];
	temp[((nx+2)*iy)+(nx+1)] = temp[((nx+2)*iy)+nx];
    }

    for(ix=0;ix<=nx+1;++ix){
	temp[ix] = temp[(nx+2)+ix];
	temp[((nx+2)*(ny+1))+ix] = temp[((nx+2)*ny)+ix];
    }
}

/*
 * Perform time evolution using FTCS FD method
 */

void evolve(int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT *temp, MYFLOAT alpha){
    
    MYFLOAT dx, dy;
    int ix, iy, ix_start, ix_end, iy_tmp;

    int th_Id = 0, numThreads = 1, block_size;

    dx = lx/nx;
    dy = ly/ny;

#pragma omp parallel private( iy_tmp, iy, ix, th_Id, ix_start, ix_end )
    {
      
#ifdef __OPENMP
      th_Id = omp_get_thread_num();
      
#pragma omp single
      {
	numThreads = omp_get_num_threads();
	block_size = nx/numThreads; 
      }
#endif
      
      ix_start =  block_size * th_Id + 1;
      ix_end = ix_start + block_size - 1;

      for( iy_tmp = 1; iy_tmp <= ny + numThreads - 1; ++iy_tmp){

	iy = iy_tmp - th_Id;

	if( ( iy >= 1 ) && ( iy <= ny ) ){

	  for(ix=ix_start;ix<=ix_end;++ix){

	    temp[((nx+2)*iy)+ix] = temp[((nx+2)*iy)+ix] + alpha*dt*
	      ( (temp[((nx+2)*(iy+1))+ix] + temp[((nx+2)*(iy-1))+ix] -2.0* temp[((nx+2)*iy)+ix])/(dy*dy) + 
		(temp[((nx+2)*iy)+(ix+1)] + temp[((nx+2)*iy)+(ix-1)] -2.0* temp[((nx+2)*iy)+ix])/(dx*dx) );
	  }
	}
#pragma omp barrier
      }    
    }
}


int main(int argc, char* argv[]){


    int i, nx, ny, n_steps, nRows, nCols, frame;
    MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
    MYFLOAT *temp;
    MYFLOAT norm, norm_ini, bound;
    FILE *fp;
    double t1, t2;

    if (argc !=3) 
       {
       printf(" FTCS finite difference solution of the heat equation \n\n");
       printf(" Usage : %s n_steps dt \n", argv[0]);
       return 1;
       }

    // number of points in the x directions
    nx=8192;
    nCols= nx + 2;
    // number of points in the y directions
    ny=8192;
    nRows= ny + 2;
    // size of the system in the x direction
    lx=2.0;
    // size of the system in the y direction
    ly=2.0;
    //center of the initial gaussian distribution    
    x0=0.0;
    y0=0.0;
    // amplitude of the initial gaussian distribution 
    sigmax=0.5;
    sigmay=0.2;
    // alpha coefficient (diffusivity) of the heat equation 
    alpha=1.0;
    // number of time steps and amplitude thereof. Command line argumentes 
    n_steps=atoi(argv[1]);
    dt=atof(argv[2]);

    temp = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));
    
    // fill temperaature array with initial condition, imposing flat boundary conditions
    init(temp, nx, ny, lx, ly, x0, y0, sigmax, sigmay);
    update_boundaries_FLAT(nx, ny,temp);
    
    printf(" Starting time evolution... \n\n ");
    t1 = cclock();
    for(i=1; i<=n_steps; ++i) {

      evolve(nx, ny, lx, ly, dt, temp, alpha);
      update_boundaries_FLAT(nx, ny, temp);

    }
    t2 = cclock();
    
    fprintf( stdout, "\n Time to solution: %.3g\n\n", t2 - t1 );

    return 0;
}

