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
#include <mpi.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#ifdef LINK_OMP
#include <omp.h>
#endif

#define PI 3.14159265359

#ifdef __SINGLE

typedef float MYFLOAT;
#define MY_MPI_FLOAT MPI_FLOAT

#else

typedef double MYFLOAT;
#define MY_MPI_FLOAT MPI_DOUBLE

#endif

#define DUMP 1000

/*
 * conversions from discrete to real coordinates
 */
MYFLOAT ix2x(int ix, int nx, MYFLOAT lx){
    return ((ix-1) - nx / 2.0)*lx/nx;
}

MYFLOAT iy2y(int iy, int ny, MYFLOAT ly){ 
    return ((iy-1) - ny / 2.0)*ly/ny;
}

/* 
 * Integral of the temperature field 
 */
MYFLOAT integral(int nx, int ny, int ny_global, MYFLOAT* val, MYFLOAT lx, MYFLOAT ly){
	MYFLOAT sum=0.0;
	int ix,iy;
    	for(iy=1;iy<=ny;++iy)
		for(ix=1;ix<=nx;++ix){
          		sum+=val[((nx+2)*iy)+ix];
 		} 
 	return(sum*lx/nx*ly/ny_global);
}


/* 
 * initialize the system with a gaussian temperature distribution. 
 * The center of the gaussian (x0,y0) as well as the amplitudes sigmax, sigmay are input parameters 
 */
void init(MYFLOAT *temp, int MyID, int nx, int ny, int ny_global, int offset, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay){ 
  int ix, iy, iy_global;
  MYFLOAT x, y;
#ifdef LINK_OMP
#pragma omp parallel for private(iy_global, ix, x, y)
#endif
  for(iy=1;iy<=ny;++iy){
    iy_global = iy + offset + ny * MyID;
    for(ix=1;ix<=nx;++ix){
      x=ix2x(ix, nx, lx);
      y=iy2y(iy_global, ny_global, ly);
      temp[((nx+2)*iy)+ix] = 1.0/(2.*PI*sigmax*sigmay)*exp(-(x-x0)*(x-x0)/(2.0*(sigmax*sigmax)) - (y-y0)*(y-y0)/(2.0*(sigmay*sigmay)) );
    }
  }
  if(MyID == 0)
    printf("--  Max (for gnuplot) -- %f\n", 1.0/(2.*PI*sigmax*sigmay) );
}

/*
 * save the temperature distribution
 * the ascii format is suitable for splot gnuplot function
 */
void save_gnuplot(FILE* fp, MYFLOAT *temp, int nx, int ny, int ny_global, MYFLOAT lx, MYFLOAT ly, int MyID, int NPE, int rest) {
    
  int ix, iy, process;
  int tag = 42;

  if(MyID == 0){
    int iy_global, offset = 0;
    MYFLOAT *tmp_buf;
    

    tmp_buf = (MYFLOAT*)malloc((nx + 2) * (ny + 2) * sizeof(MYFLOAT));

    // the double newline will be interpreted by gnuplot as a new data block
    fprintf(fp, "\n\n");
    for(iy = 1; iy <= ny; ++iy){
      iy_global = iy + ny * MyID;
      for(ix = 1; ix <= nx; ++ix)
	fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy_global,ny_global,ly), temp[((nx+2)*iy)+ix]);
    }

    /* recive data from other processes */
    for(process = 1; process < NPE; process++){
      if(rest != 0 && process == rest){
	ny--;
	offset = rest;
      }
      MPI_Recv(tmp_buf, (nx + 2) * (ny + 2), MY_MPI_FLOAT, process, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(iy = 1; iy <= ny; ++iy){
	iy_global = iy + offset + ny * process;
	for(ix = 1; ix <= nx; ++ix)
	  fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy_global,ny_global,ly), tmp_buf[((nx+2)*iy)+ix]);
      }
    }
    if(rest != 0)
      ny++;

    free(tmp_buf);
  }
  else{
    MPI_Send(temp, (nx + 2) * (ny + 2), MY_MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
  }
}

/*
 * Updating boundaries. Boundaries are flat: copying the value of the neighbour
 */

void update_boundaries_FLAT(int MyID, int NPE, int nx, int ny, MYFLOAT *temp){

  /* int ix, iy,  */
  int prev, next;

  int send_up_tag = 42;
  int send_down_tag = 24;

  MPI_Request send_up, send_down, rec_up;
  MPI_Status status, rec_down;

  prev = (MyID + NPE - 1) % NPE;
  next = (MyID + 1) % NPE;

  if(MyID == 0)
    prev = MPI_PROC_NULL;
  if(MyID == NPE - 1)
    next = MPI_PROC_NULL;

  MPI_Isend(temp + nx + 2, (nx+2), MY_MPI_FLOAT, prev, send_up_tag, MPI_COMM_WORLD, &send_up);
  MPI_Isend(temp + (nx + 2) * ny, (nx+2), MY_MPI_FLOAT, next, send_down_tag, MPI_COMM_WORLD, &send_down);
    
  MPI_Irecv(temp, (nx+2), MY_MPI_FLOAT, prev, send_down_tag, MPI_COMM_WORLD, &rec_up);
  MPI_Recv(temp + (nx + 2) * (ny + 1), (nx+2), MY_MPI_FLOAT, next, send_up_tag, MPI_COMM_WORLD, &rec_down);

  MPI_Wait(&rec_up, &status);  
  MPI_Wait(&send_up, &status);
  MPI_Wait(&send_down, &status);
}

/*
 * Perform time evolution using FTCS FD method
 */

void evolve(int nx, int ny, int ny_global, MYFLOAT lx, MYFLOAT ly, MYFLOAT dt, MYFLOAT *temp, MYFLOAT *temp_new, MYFLOAT alpha){
    
    MYFLOAT dx, dy;
    int ix, iy;

    dx = lx/nx;
    dy = ly/ny_global;
#ifdef LINK_OMP
#pragma omp parallel for private(ix)
#endif
    for(iy=1;iy<=ny;++iy)
	for(ix=1;ix<=nx;++ix){
	    temp_new[((nx+2)*iy)+ix] = temp[((nx+2)*iy)+ix] + alpha*dt*
              ( (temp[((nx+2)*(iy+1))+ix] + temp[((nx+2)*(iy-1))+ix] -2.0* temp[((nx+2)*iy)+ix])/(dy*dy) + 
                (temp[((nx+2)*iy)+(ix+1)] + temp[((nx+2)*iy)+(ix-1)] -2.0* temp[((nx+2)*iy)+ix])/(dx*dx) );
	}
#ifdef LINK_OMP
#pragma omp parallel for private(ix)
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


    int i, nx, ny, n_steps, nRows, nCols;
    MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
    MYFLOAT *temp, *temp_new;
    MYFLOAT norm, norm_ini, bound, norm_ini_root, norm_root;
    FILE *fp;
    
    double t_start, t_end, t_sum = 0;

    int MyID, NPE, ny_global, rest, offset = 0;

    if (argc !=3) 
       {
       printf(" FTCS finite difference solution of the heat equation \n\n");
       printf(" Usage : %s n_steps dt \n", argv[0]);
       return 1;
       }

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
    MPI_Comm_size(MPI_COMM_WORLD, &NPE);

    // number of points in the x directions
    nx = 2048*4;
    nCols= nx + 2;
    // number of points in the y directions
    ny_global = 2048*4;

    if(ny_global < NPE){
      printf("\n\tNumber of rows must be greater than NPE\n\tExit!\n");
      return 1;
    }

    ny = ny_global / NPE;
    rest = ny_global % NPE;

    if(rest != 0 && MyID < rest)
      ny++;
    else
      offset = rest;

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
    n_steps = atoi(argv[1]);
    dt = atof(argv[2]);

    // checking stability of the method 
    bound = 1.0/(2.0*alpha)*(1.0/(nx*nx/(lx*lx)+ny*ny/(ly*ly))); 
    if(dt > bound)
      if(MyID == 0)
	printf("Warning, unstable regime, reduce dt \n");

    // HINT - the exercise requires to work on distributed data
    temp = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));
    temp_new = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));

    // fill temperaature array with initial condition, imposing flat boundary conditions
    init(temp, MyID, nx, ny, ny_global, offset, lx, ly, x0, y0, sigmax, sigmay);
    update_boundaries_FLAT(MyID, NPE, nx, ny,temp);

    norm_ini=integral(nx, ny, ny_global, temp, lx, ly);
    norm_ini_root = 0.;
    norm_root = 0.;
    MPI_Reduce(&norm_ini, &norm_ini_root, 1, MY_MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(MyID == 0){
      printf("\tInitial integral value: %f\n", norm_ini_root);
      fp = fopen("heat_diffusion.dat", "w");
      printf("\tStarting time evolution... \n\n ");
    }

    if( (nx * ny ) < 10001 )
      save_gnuplot(fp, temp, nx, ny, ny_global, lx, ly, MyID, NPE, rest);
   

    t_start = seconds();
    
    for(i=1; i<=n_steps; ++i) {
      // performing TFSC-FD step and updating boundaries
      evolve(nx, ny, ny_global, lx, ly, dt, temp, temp_new, alpha);
      update_boundaries_FLAT(MyID, NPE, nx, ny, temp);
      
      //saving a snapshot every DUMP time steps
      if( i % DUMP == 0 && (nx * ny) < 10001){
	t_end = seconds();
	save_gnuplot(fp, temp, nx, ny, ny_global, lx, ly, MyID, NPE, rest);
	t_sum += t_end - t_start;
	t_start = seconds();
      }
    }
    
    t_end = seconds();
    t_sum += t_end- t_start;

    // checking final norm (~total energy)
    norm=integral(nx, ny, ny_global, temp, lx, ly);
    MPI_Reduce(&norm, &norm_root, 1, MY_MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(MyID == 0){
      printf("\tIntegral end: %f\n",norm_root);
      printf("\tTime used: %lg s\n", t_sum);
      fclose(fp);
    }
    
    MPI_Finalize();
    
    return 0;
}

