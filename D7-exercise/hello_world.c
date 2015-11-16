#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char** argv){
  if(argc > 1)
    int nt = atoi(argv[1]);
  else
    nt = 0;

  omp_set_num_threads(nt);
  printf("\n\n\tSequential part...\n");
#pragma omp parallel
{
  printf("\thello from thread number %d\n", omp_get_thread_num());

}
  printf("\n\tBack to sequential part...\n\n");

  return 0;
}
