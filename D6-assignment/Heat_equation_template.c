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

#define PI 3.14159265359

#ifdef __SINGLE

typedef float MYFLOAT;
#define MY_MPI_FLOAT MPI_FLOAT

#else

typedef double MYFLOAT;
#define MY_MPI_FLOAT MPI_DOUBLE

#endif

#define DUMP 100


/* 
 * conversions from discrete to real coordinates
 * HINT -  Consider the shift from local coordinates to global coodinates
 */
MYFLOAT ix2x(int ix, int nx, MYFLOAT lx){
    return ((ix-1) - nx / 2.0)*lx/nx;
}

MYFLOAT iy2y(int iy, int ny, MYFLOAT ly){ 
    return ((iy-1) - ny / 2.0)*ly/ny;
}


/* Integral of the temperature field 
 *  HINT - In the parallel case, the integration on each process will yeld a the integral on the local region
 */
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
 * initialize the system with a gaussian temperature distribution. 
 * The center of the gaussian (x0,y0) as well as the amplitudes sigmax, sigmay are input parameters 
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
    printf("--  Max (for gnuplot) -- %f\n", 1.0/(2.*PI*sigmax*sigmay) );
    return 0;
}

/*
 * save the temperature distribution
 * the ascii format is suitable for splot gnuplot function
 * HINT - print distributed data from ROOT_PROCESS
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
 * HINT -  handle the ghost-cell exchange 
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

void evolve(int nx, int ny, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT *temp, MYFLOAT *temp_new, MYFLOAT alpha){
    
    MYFLOAT dx, dy;
    int ix, iy;

    dx = lx/nx;
    dy = ly/ny;

    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
	    temp_new[((nx+2)*iy)+ix] = temp[((nx+2)*iy)+ix] + alpha*dt*
              ( (temp[((nx+2)*(iy+1))+ix] + temp[((nx+2)*(iy-1))+ix] -2.0* temp[((nx+2)*iy)+ix])/(dy*dy) + 
                (temp[((nx+2)*iy)+(ix+1)] + temp[((nx+2)*iy)+(ix-1)] -2.0* temp[((nx+2)*iy)+ix])/(dx*dx) );
	}


    for(iy=0;iy<=ny+1;++iy)
	for(ix=0;ix<=nx+1;++ix)
	    temp[((nx+2)*iy)+ix] = temp_new[((nx+2)*iy)+ix];

}



int main(int argc, char* argv[]){


    int i, nx, ny, n_steps, nRows, nCols;
    MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
    MYFLOAT *temp, *temp_new;
    MYFLOAT norm, norm_ini, bound;
    FILE *fp;

    if (argc !=3) 
       {
       printf(" FTCS finite difference solution of the heat equation \n\n");
       printf(" Usage : %s n_steps dt \n", argv[0]);
       return 1;
       }

    // number of points in the x directions
    nx=100;
    nCols= nx + 2;
    // number of points in the y directions
    ny=50;
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

    // checking stability of the method 
    bound = 1.0/(2.0*alpha)*(1.0/(nx*nx/(lx*lx)+ny*ny/(ly*ly))); 
    if (dt>bound  )
       {
       printf("Warning, unstable regime, reduce dt \n");
       }

    // HINT - the exercise requires to work on distributed data
    temp = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));
    temp_new = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));

    // fill temperaature array with initial condition, imposing flat boundary conditions
    init(temp, nx, ny, lx, ly, x0, y0, sigmax, sigmay);
    update_boundaries_FLAT(nx, ny,temp);

    // calculating initial norm (~ total energy)
    // HINT - perform a parallel integration
    norm_ini=integral(nx, ny, temp, lx, ly);
    printf(" Initial integral value: %f\n", norm_ini);

    fp = fopen("heat_diffusion.dat", "w");

    printf(" Starting time evolution... \n\n ");

    if( (nx * ny ) < 10001 )
      save_gnuplot(fp, temp, nx, ny, nly, lx, ly, rank, nproc);

    for(i=1; i<=n_steps; ++i) {

         // performing TFSC-FD step and updating boundaries
         evolve(nx, ny, lx, ly, dt, temp, temp_new, alpha);
         update_boundaries_FLAT(nx, ny, temp);

         // saving a snapshot every DUMP time steps
	 if( i % DUMP == 0 && (nx * ny ) < 10001 )
	   save_gnuplot(fp, temp, nx, ny, nly, lx, ly, rank, nproc); 

    }

    // checking final norm (~total energy)
    norm=integral(nx, ny, temp, lx, ly); 
    printf(" Integral end: %f\n",norm);

    fclose(fp); 

    return 0;
}

