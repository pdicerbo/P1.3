/*
 *     Using Time Forward Space Centered finite difference method this code solves the heat equation
 *
 *                            u,t = alpha* (u,xx + u,yy) 
 *
 *    using a domain decomposition and MPI communications.
 *    Boundaries are flat: copying the value of the neighbour
 *    Created by G.P. Brandino for MHPC
 *    Last Revision: October 2015
 */

#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265359

#ifdef __SINGLE

typedef float MYFLOAT;
#define MY_MPI_FLOAT MPI_FLOAT

#else

typedef double MYFLOAT;
#define MY_MPI_FLOAT MPI_DOUBLE

#endif

#define DUMP 100


double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

void swap_buffer( MYFLOAT ** send, MYFLOAT ** recv){

  MYFLOAT * tmp;

  tmp = (* send);
  (* send) = (* recv);
  (* recv) = tmp;
}


/* 
 * conversions from discrete to real coordinates
 */
MYFLOAT ix2x(int ix, int nx, MYFLOAT lx){
    return ((ix-1) - nx / 2.0)*lx/nx;
}

MYFLOAT iy2y(int iy, int ny, MYFLOAT ly, int offset){
 
    return ((iy-1) - ny / 2.0 + offset )*ly/ny;
}


/* Integral of temperature field */
MYFLOAT integral(int nx, int ny, int nly, MYFLOAT* val, MYFLOAT lx, MYFLOAT ly){
	MYFLOAT sum=0.0;
	int ix,iy;
    	for(iy=1;iy<=nly;++iy)
		for(ix=1;ix<=nx;++ix){
          		sum+=val[((nx+2)*iy)+ix];
 		} 
 	return(sum*lx/nx*ly/ny);
}


/*
 * initialize the system with a gaussian temperature distribution. The center of the gaussian (x0,y0) as well as the amplitudes sigmax, sigmay are input parameters 
 */
void init(MYFLOAT *temp, int nx, int ny, int nly, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay, int rank, int nproc){ 
    int i,ix,iy, nly_r;
    int offset;
    MYFLOAT x,y;

    if (rank < ny%nproc)
	offset = (ny/nproc+1)*rank;
    else 
	offset = ny/nproc*rank + ny%nproc;

    for(iy=1;iy<=nly;++iy)
	for(ix=1;ix<=nx;++ix){
	    x=ix2x(ix, nx, lx);
	    y=iy2y(iy, ny, ly, offset);
	    temp[((nx+2)*iy)+ix] = 1.0/(2.*PI*sigmax*sigmay)*exp(-(x-x0)*(x-x0)/(2.0*(sigmax*sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay*sigmay)) );
	}
    if (rank ==0 )
	printf("--  Max (for gnuplot) -- %f\n", 1.0/(2.*PI*sigmax*sigmay) );
}

/*
 * save the temperature distribution
 * the ascii format is suitable for splot gnuplot function
 */
int save_gnuplot(FILE* fp, MYFLOAT *temp, int nx, int ny, int nly, MYFLOAT lx, MYFLOAT ly, int rank, int nproc) {
    
    int i,ix,iy, proc;
    int nly_p[nproc], offset, nly_r;
    MYFLOAT* buf;
    MPI_Status status;

    for (i=0; i<nproc;i++){
            nly_p[i] = ny/nproc;
            if (i < (ny % nproc) )
                nly_p[i]++;
    }

    if (rank == 0){
	buf=(MYFLOAT*)malloc((nx+2)*(ny+2)*sizeof(MYFLOAT));
    	fprintf(fp, "\n\n");
    	for(iy=1;iy<=nly;++iy){		
		for(ix=1;ix<=nx;++ix)
	    	   fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy,ny,ly,0), temp[((nx+2)*iy)+ix]);
    	}
        for(proc=1; proc<nproc; proc++){

		if (proc < ny%nproc)
        		offset = (ny/nproc+1)*proc;
    		else
       			offset = ny/nproc*proc + ny%nproc;
	
		MPI_Recv(buf, (nx+2)*(nly_p[proc]+2), MY_MPI_FLOAT, proc, 0, MPI_COMM_WORLD, &status);
                for(iy=1;iy<=nly_p[proc];++iy){
              		  for(ix=1;ix<=nx;++ix)
	                	fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy,ny,ly,offset), buf[((nx+2)*iy)+ix]);   
			  }
        }
        free(buf);
    }
    else
    {
    MPI_Send(temp, (nx+2)*(nly+2), MY_MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }
    
    return 0;
}

void update_boundaries_FLAT(int nx, int nly, MYFLOAT *temp, int rank, int nproc, int rank_up, int rank_down){

    int ix, iy;
    MPI_Status status, status1;


    if ( rank == 0 ){
           for(ix=0;ix<=nx+1;++ix)
               temp[ix] = temp[(nx+2)+ix];
    }

    if ( rank == nproc-1 ){
           for(ix=0;ix<=nx+1;++ix)
		temp[((nx+2)*(nly+1))+ix] = temp[((nx+2)*nly)+ix];
    }


    MPI_Sendrecv(&temp[(nx+2)], nx+2, MY_MPI_FLOAT, rank_down, 0, &temp[((nx+2)*(nly+1))], nx+2, MY_MPI_FLOAT, rank_up, 0, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&temp[((nx+2)*(nly))], nx+2, MY_MPI_FLOAT, rank_up, 0, &temp[0], nx+2, MY_MPI_FLOAT, rank_down, 0, MPI_COMM_WORLD, &status1);

    for(iy=0;iy<=nly+1;++iy){
        temp[(nx+2)*iy] = temp[((nx+2)*iy)+1];
        temp[((nx+2)*iy)+(nx+1)] = temp[((nx+2)*iy)+nx];
    }


}


void evolve(int nx, int ny, int nly, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT **temp_p, MYFLOAT **temp_new_p, MYFLOAT alpha){
    
    MYFLOAT dx, dy;
    int ix, iy;
    MYFLOAT *temp, *temp_new;

    dx = lx/nx;
    dy = ly/ny;

    temp=*temp_p;
    temp_new=*temp_new_p;

    for(iy=1;iy<=nly;++iy)
	for(ix=1;ix<=nx;++ix){
	    temp_new[((nx+2)*iy)+ix] = temp[((nx+2)*iy)+ix] + alpha*dt*
              ( (temp[((nx+2)*(iy+1))+ix] + temp[((nx+2)*(iy-1))+ix] -2.0* temp[((nx+2)*iy)+ix])/(dy*dy) + 
                (temp[((nx+2)*iy)+(ix+1)] + temp[((nx+2)*iy)+(ix-1)] -2.0* temp[((nx+2)*iy)+ix])/(dx*dx) );
	}

    swap_buffer( temp_p, temp_new_p );
}



int main(int argc, char* argv[]){


    int i, nx, ny, nly, n_steps, nRows, nCols;
    MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
    MYFLOAT *temp, *temp_new;
    MYFLOAT norm, norm_l, norm_ini, norm_ini_l, bound;
    FILE *fp;
    int rank, nproc, rank_up, rank_down;

    double t_start, t_stop;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (argc !=3) {
      if ( rank==0){
	
	printf(" FTCS finite difference solution of the heat equation \n\n");
	printf(" Usage : %s n_steps dt \n", argv[0]);
      }
      MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    
    // number of points in the x directions
    nx=4000;
    // number of points in the y directions
    ny=4000;
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
    if (dt>bound && rank ==0  )
       {
       printf("Warning, unstable regime, reduce dt \n");
       }

    if (ny < nproc ) {
         if (rank==0) {
            printf("\nTrying to run with %d processes, while the maximum allowed size is %d (NY)\n",nproc,ny);
         }
         MPI_Finalize();
    }

    // balacing workload
    nly = ny/nproc;
    if (rank < (ny % nproc) )
         nly++;

    nCols= nx+2;
    nRows= nly+2;


    temp = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));
    temp_new = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));

    rank_down = rank-1;
    rank_up = rank+1;

    if(rank_down < 0) rank_down = MPI_PROC_NULL;
    if(rank_up == nproc) rank_up = MPI_PROC_NULL;


    init(temp, nx, ny, nly, lx, ly, 0.0, 0.0, sigmax, sigmay, rank, nproc);
    update_boundaries_FLAT(nx, nly,temp, rank, nproc, rank_up, rank_down);

    norm_ini_l=integral(nx, ny, nly, temp, lx, ly);
    norm_ini=0.0;
    MPI_Reduce(&norm_ini_l,&norm_ini,1, MY_MPI_FLOAT,MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank ==0){  
       printf(" Initial integral value: %f\n", norm_ini);
       fp = fopen("heat_diffusion.dat", "w");
       printf(" Starting time evolution... \n\n ");
    }

    if( (nx * ny ) < 10001 )
      save_gnuplot(fp, temp, nx, ny, nly, lx, ly, rank, nproc);
    
    t_start = seconds();
    for(i=1; i<=n_steps; ++i) {

         evolve(nx, ny, nly, lx, ly, dt, &temp, &temp_new, alpha);
         update_boundaries_FLAT(nx, nly,temp, rank, nproc, rank_up, rank_down); 
	 if( i % DUMP == 0 && (nx * ny ) < 10001 )
	   save_gnuplot(fp, temp, nx, ny, nly, lx, ly, rank, nproc); 
    }
    t_stop = seconds();

    norm_l=integral(nx, ny, nly, temp, lx, ly); 
    norm=0.0;
    MPI_Reduce(&norm_l,&norm,1, MY_MPI_FLOAT,MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank==0){
	    printf(" Integral end: %f\n",norm);
	    fprintf( stdout, "\nTime to solution: %.3g", t_stop - t_start );
	    fclose(fp); 
    }

    MPI_Finalize();
    return 0;
}

