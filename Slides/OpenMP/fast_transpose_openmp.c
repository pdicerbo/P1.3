/*
 *    Copyright (C) 2001-2014 The Abdus Salam, International Centre of Theoretical Physics (ICTP)
 *
 *    This file is distributed under the terms of the GNU General Public License. 
 *    See http://www.gnu.org/copyleft/gpl.txt 
 *
 *    This code implements a multi-threaded version of the fast-transpose problem presented in class
 *    during the course of "Computer Archietctures for HPC".
 *    The code is parallelized using OpenMP and it is aimed for teaching purpose. 
 *    Parallelization is applied for initializzation and to transpose the matrix.
 *    Transposition is performed out of place (input and output stored on different buffers).
 *    The program works only for square Matrixes of dimention multiple of the Block dimention.
 *    The cache block is allocated as a statuc array buffer of BLOCK_SIZE dimension
 *
 *    Created by I. Girotto for MHPC
 *    Last Revision: Novemeber 2015
 */


#include <stdlib.h>
#include <stdio.h>

#ifdef __OPENMP
#include <omp.h>
#endif

#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>

#define BLOCK_SIZE 64 // store the dimension of the cache block 

double cclock()
/* Returns elepsed seconds past from the last call to timer rest */
{

  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

int main( int argc, char * argv [] ) {

  double *A,  *AT;
  double t1, t2;
  int im, jm, ib, jb, ioff, joff, num_blocks;
  int m_size, b_size, nthreads;
  double buffer[ BLOCK_SIZE * BLOCK_SIZE];

  if( argc < 2 ){
    fprintf( stderr, "Error. The program runs as following: %s [MATRIXDIM].\nProgram exit ...\n", argv[0] );
    exit(1);
  }
  
  m_size = atoi(argv[1]); // store the dimension of the Matrices A & AT 

  //compute the number of block
  num_blocks = m_size / BLOCK_SIZE;

  if( m_size % BLOCK_SIZE ){
    fprintf( stderr, "Error. Inconsistent parameters. BLOCKSIZE = %d is not factor of MATRIXDIM\nProgram exit ...\n", BLOCK_SIZE);
    exit(1);
  }


  // Allocation of the  Matrix memory
  A = (double*) malloc( m_size * m_size * sizeof(double) );
  AT = (double*) malloc( m_size * m_size * sizeof(double) );

  //Initialization of the A matrix
#ifdef __OPENMP
#pragma omp parallel for private(im, jm)
#endif  
  for( im = 0; im < m_size; im++ )
    for( jm = 0; jm < m_size; jm++ ){
      A[ im * m_size + jm ] = (double) ( ( im * m_size ) + jm );
    }
  
  
  t1 = cclock();
  
  /* Implement loops over the blocks of the main Matrices while transposing internally to the cache clock */
#ifdef __OPENMP
#pragma omp parallel private(im, jm, ib, jb, ioff, joff, buffer) firstprivate(num_blocks)
{

# pragma omp for collapse(2)
#endif
  for( ioff = 0; ioff < m_size; ioff += BLOCK_SIZE )
    {
      for( joff = 0; joff < m_size; joff += BLOCK_SIZE)
	{
	  // copy from A to buffer block
	  for( ib = 0, im = ioff; ib < BLOCK_SIZE; ib++, im++)
	    for (jb = 0, jm = joff; jb < BLOCK_SIZE; jb++, jm++)
	      buffer[ ib * BLOCK_SIZE + jb ] = A[ im * m_size + jm];
	  
	  // copy from a transposed buffer to AT	
	  for( ib = 0, im = joff; ib < BLOCK_SIZE; ib++, im++)
	    for (jb = 0, jm = ioff; jb < BLOCK_SIZE; jb++, jm++)
	      AT[ im * m_size + jm ] = buffer [ jb * BLOCK_SIZE + ib];
	  
	} 
    }
#ifdef __OPENMP
 } /* omp parallel */
#endif
 
 t2 = cclock();
 
#ifdef __OPENMP
#pragma omp parallel
 { 
#pragma omp single
    fprintf( stdout, "\n\tOpenMP is enabled. Run executed with %d threads.", omp_get_num_threads() );
 }
#endif
 fprintf( stdout, "\n\tTime to solution: %.3g\n\n", t2 - t1 );
 
 free(A);
 free(AT);
 
 return 0;
}
