/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%
   %%%% Implement a recursive algorithm for matrix transposition
   %%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

void transpose(int* A, int n, int m, MPI_Comm, comm);

int main(int argc,char **argv) {

  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int nprocs, procno;
  MPI_Comm_rank(comm,&procno);
  MPI_Comm_size(comm,&nprocs);

  //defining the size of the matrix and initializing
  int n=2
  int m=2
    
  int* A = (int*)malloc(n * m * sizeof(int));
  
  for (int i=0, i<n, i++) {
    for (int j=0, j<m, j++) {
      A[i * m + j] = i * m + j;
    }
  }

  transpose(A, n, m, MPI_COMM_WORLD);
    
  if (procno == 0) {
     printf("Transposed matrix:\n");
        for (int i = 0; i < m; i++) {
           for (int j = 0; j < n; j++) {
              printf("%d ", A[j * m + i]);
            }
            printf("\n"); 
        }
  }
   
    MPI_Finalize();
    return 0;
}

void transpose(int* A, int n, int m, MPI_Comm, comm) {
  
  int nprocs, procno;
  MPI_Comm_rank(comm,&procno);
  MPI_Comm_size(comm,&nprocs);
  
  //incomplete at the moment
  
  
  
  
  
  
  
  
