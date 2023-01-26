/****************************************************************
 ****
 **** This program file is part of the book 
 **** `Parallel programming for Science and Engineering'
 **** by Victor Eijkhout, eijkhout@tacc.utexas.edu
 ****
 **** copyright Victor Eijkhout 2012-8
 ****
 **** MPI Exercise for prime number factorization. So to speak.
 ****
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc,char **argv) {

  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int nprocs, procno;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&procno);

  int bignum = 2000000111, maxfactor = 45200;

  int range = maxfactor / nprocs;
  int start = range * procno + 2;
  int end = start + range - 1;

  for (int myfactor = start; myfactor <= end; myfactor++) {
    if (bignum % myfactor == 0) {
      printf("Processor %d found factor %d\n", procno, myfactor);
    }
  }

  MPI_Finalize();
  return 0;
}          
