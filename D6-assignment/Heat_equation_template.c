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
void init(MYFLOAT *temp, int MyID, int nx, int ny, int ny_global, int offset, MYFLOAT lx, MYFLOAT ly, MYFLOAT x0, MYFLOAT y0, MYFLOAT sigmax, MYFLOAT sigmay){ 
  int ix, iy, iy_global;
    MYFLOAT x, y;

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
 * HINT - print distributed data from ROOT_PROCESS
 */
void save_gnuplot(FILE* fp, MYFLOAT *temp, int nx, int ny, int ny_global, MYFLOAT lx, MYFLOAT ly, int MyID, int NPE, int rest) {
    
  int ix, iy, process;
  int tag = 42;

  if(MyID == 0){
    MYFLOAT *tmp_buf;
    tmp_buf = (MYFLOAT*)malloc((nx + 2) * (ny + 2) * sizeof(MYFLOAT));

    // the double newline will be interpreted by gnuplot as a new data block
    fprintf(fp, "\n\n");
    for(iy = 1; iy <= ny; ++iy){
      for(ix = 1; ix <= nx; ++ix)
	fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy,ny,ly), temp[((nx+2)*iy)+ix]);
    }
    for(process = 1; process < NPE; process++){
      if(rest != 0 && process == rest)
	ny--;
      MPI_Recv(tmp_buf, (nx + 2) * (ny + 2), MY_MPI_FLOAT, process, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(iy = 1; iy <= ny; ++iy){ 
	for(ix = 1; ix <= nx; ++ix)
	  fprintf(fp, "\t%f\t%f\t%g\n", ix2x(ix,nx,lx), iy2y(iy,ny,ly), tmp_buf[((nx+2)*iy)+ix]);
      }
    }
    if(rest != 0)
      ny++;

  }
  else{
    MPI_Send(temp, (nx + 2) * (ny + 2), MY_MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
  }
}

/*
 * Updating boundaries. Boundaries are flat: copying the value of the neighbour
 * HINT -  handle the ghost-cell exchange 
 */

void update_boundaries_FLAT(int MyID, int NPE, int nx, int ny, MYFLOAT *temp){

  int ix, iy, prev, next;
  MYFLOAT *buf_send_up, *buf_send_down, *buf_rec_up, *buf_rec_down;

  int send_up_tag = 42;
  int send_down_tag = 24;

  MPI_Request send_up, send_down, rec_up;
  MPI_Status status, rec_down;

  prev = (MyID + NPE - 1) % NPE;
  next = (MyID + 1) % NPE;

  for(iy = 1; iy <= ny; iy++){
    temp[(nx + 2)*iy] = temp[(nx + 2)*iy + 1];
    temp[(nx + 1) + (nx + 2)*iy] = temp[nx + (nx + 2)*iy];
  }

  if(MyID > 0 && MyID < NPE - 1){

    buf_send_up   = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    buf_send_down = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    buf_rec_up    = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    buf_rec_down  = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
      
    for(ix = 0; ix < nx + 2; ix++){
      buf_send_up[ix] = temp[nx + 2 + ix];
      buf_send_down[ix] = temp[(nx + 2) * ny + ix];
    }

    MPI_Isend(buf_send_up, (nx+2), MY_MPI_FLOAT, prev, send_up_tag, MPI_COMM_WORLD, &send_up);
    MPI_Isend(buf_send_down, (nx+2), MY_MPI_FLOAT, next, send_down_tag, MPI_COMM_WORLD, &send_down);
    
    MPI_Irecv(buf_rec_up, (nx+2), MY_MPI_FLOAT, prev, send_down_tag, MPI_COMM_WORLD, &rec_up);
    MPI_Recv(buf_rec_down, (nx+2), MY_MPI_FLOAT, next, send_up_tag, MPI_COMM_WORLD, &rec_down);

    for(ix = 0; ix < nx + 2; ix++)
      temp[(nx + 2) * (ny + 1) + ix] = buf_rec_down[ix];

    MPI_Wait(&rec_up, &status);
    
    for(ix = 0; ix < nx + 2; ix++)
      temp[ix] = buf_rec_up[ix];

    MPI_Wait(&send_up, &status);
    MPI_Wait(&send_down, &status);

    free(buf_send_up);
    free(buf_send_down);
    free(buf_rec_up);
    free(buf_rec_down);
  }
  else if(MyID == 0){

    buf_send_down = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    buf_rec_down  = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
      
    for(ix = 0; ix < nx + 2; ix++)
      buf_send_down[ix] = temp[(nx + 2) * ny + ix];

    MPI_Isend(buf_send_down, (nx+2), MY_MPI_FLOAT, next, send_down_tag, MPI_COMM_WORLD, &send_down);

    MPI_Irecv(buf_rec_down, (nx+2), MY_MPI_FLOAT, next, send_up_tag, MPI_COMM_WORLD, &rec_up);

    /* update first row of matrix */
    for(iy = 0; iy < nx + 2; iy++)
      temp[iy] = temp[(nx + 2) + iy];

    MPI_Wait(&rec_up, &status);

    for(ix = 0; ix < nx + 2; ix++)
      temp[(nx + 2) * (ny + 1) + ix] = buf_rec_down[ix];

    MPI_Wait(&send_down, &status);

    free(buf_send_down);
    free(buf_rec_down);
  }
  else{
    buf_send_up = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    buf_rec_up  = (MYFLOAT*)malloc((nx + 2) * sizeof(MYFLOAT));
    
    for(ix = 0; ix < nx + 2; ix++)
      buf_send_up[ix] = temp[nx + 2 + ix];

    MPI_Isend(buf_send_up, (nx+2), MY_MPI_FLOAT, prev, send_up_tag, MPI_COMM_WORLD, &send_up);
    
    MPI_Irecv(buf_rec_up, (nx+2), MY_MPI_FLOAT, prev, send_down_tag, MPI_COMM_WORLD, &rec_up);
    
    /* update last row of matrix */
    for(ix = 0; ix < nx + 2; ix++)
      temp[(nx + 2) * (ny + 1) + ix] = temp[(nx + 2) * ny + ix];

    MPI_Wait(&rec_up, &status);
    
    for(ix = 0; ix < nx + 2; ix++)
      temp[ix] = buf_rec_up[ix];

    MPI_Wait(&send_up, &status);

    free(buf_send_up);
    free(buf_rec_up);
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

void print_lines(double* buffer, int cols, int rows){
  int i, j, n_i;
  for(i = 1; i <= rows; i++){
    n_i = i * (cols + 2);
    for(j = 1; j <= cols; j++)
      printf("\t%lg", buffer[j + n_i]);
    printf("\n");
  }
}


void std_out_print(double* tmp_buf, int MyRank, int SIZE, int block, int rest, int NPE){

  int process;

  int tag = 42;
  MPI_Status status;

  if(MyRank == 0){
    /* print matrix */
    printf("\n=============================\n");
    printf("\tPRINTING MATRIX (NPE = %d, rest = %d):\n\n", NPE, rest);
    
    print_lines(tmp_buf, SIZE, block);
    
    for(process = 1; process < NPE; process++){
      if(rest != 0 && process == rest)
	block--;
      MPI_Recv(tmp_buf, (SIZE + 2) * (block + 2), MPI_DOUBLE, process, tag, MPI_COMM_WORLD, &status);
      print_lines(tmp_buf, SIZE, block);
    }
    printf("\n=============================\n");
  }
  else{ /* MyRank > 0 */
    MPI_Send(tmp_buf, (SIZE + 2) * (block + 2), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  }
}

int main(int argc, char* argv[]){


    int i, nx, ny, n_steps, nRows, nCols;
    MYFLOAT lx, ly, alpha, dt, x0, y0, sigmax, sigmay;
    MYFLOAT *temp, *temp_new;
    MYFLOAT norm, norm_ini, bound;
    FILE *fp;

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
    nx = 100;
    nCols= nx + 2;
    // number of points in the y directions
    ny_global = 50;

    if(ny_global < NPE != 0){
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
       printf("Warning, unstable regime, reduce dt \n");

    // HINT - the exercise requires to work on distributed data
    temp = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));
    temp_new = (MYFLOAT *) malloc (nRows*nCols*sizeof(MYFLOAT));

    // fill temperaature array with initial condition, imposing flat boundary conditions
    init(temp, MyID, nx, ny, ny_global, offset, lx, ly, x0, y0, sigmax, sigmay);

    /* initialization check */
    /* std_out_print(temp, MyID, nx, ny, rest, NPE); */
    /* MPI_Finalize(); */
    /* return 0; */

    update_boundaries_FLAT(MyID, NPE, nx, ny,temp);

    // calculating initial norm (~ total energy)
    // HINT - perform a parallel integration
    norm_ini=integral(nx, ny, temp, lx, ly);
    printf(" Initial integral value: %f\n", norm_ini);

    fp = fopen("heat_diffusion.dat", "w");

    printf(" Starting time evolution... \n\n ");

    if( (nx * ny ) < 10001 )
      save_gnuplot(fp, temp, nx, ny, ny_global, lx, ly, MyID, NPE, rest);

    for(i=1; i<=n_steps; ++i) {

         // performing TFSC-FD step and updating boundaries
         evolve(nx, ny, lx, ly, dt, temp, temp_new, alpha);
         update_boundaries_FLAT(MyID, NPE, nx, ny, temp);

         //saving a snapshot every DUMP time steps
	 if( i % DUMP == 0 && (nx * ny) < 10001)
	   save_gnuplot(fp, temp, nx, ny, ny_global, lx, ly, MyID, NPE, rest);
    }

    // checking final norm (~total energy)
    norm=integral(nx, ny, temp, lx, ly); 
    printf(" Integral end: %f\n",norm);

    fclose(fp); 

    return 0;
}

