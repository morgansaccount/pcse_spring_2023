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

int main(int argc,char **argv) {

void transpose(double* A, int n, int i, int j, MPI_Comm comm);

  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  //transpose the matrix in place
  if (size == 1) {
  for (int i=0, i<n, i++) {
    for (int j=i+1, j<n, j++) {
      double temp = A[i*n+j];
      A[i*n+j] = A[j*n+i];
      A[j*n+i] = temp;
      }
    }
  return;
  }

//split the communicator into 4 parts
int comnu = rank % 4;
MPI_Comm subcomm;
MPI_Comm_split(comm, comnu, rank, $subcomm);

//divide the matric into 4 sub-matrices
int submat = n / 2;
double* A01 = A;
double* A02 = A + submat;
double* A03 = A + submat*n;
double* A04 = A + submat*n + submat;

//transpose each recursively
transpose(submat, A01, subcomm);
transpose(submat, A02, subcomm);
transpose(submat, A03, subcomm);
transpose(submat, A04, subcomm);

//combine sub-matrices
double* B = (double*)malloc(n*n * sizeof(double));
for (int i = 0; i < submat; i++) {
   for (int j = 0; j < submat; j++) {
      B[i*n+j] = A01[i*n+j];
      B[i*n+j+submat] = A02[i*n+j];
      B[(i+submat)*n+j] = A03[i*n+j];
      B[(i+submat)*n+j+submat] = A04[i*n+j];
   }
}

//copy the result
memcpy(A,B,n*n * sizeof(double));
free(B);
MPI_Comm_free)&subcomm);
}




