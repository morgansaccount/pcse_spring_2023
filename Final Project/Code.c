//Part 1: Vectorize a Matrix Multiplication Algorithm

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>

#define N 1000

void matrix_multiply(float *A, float *B, float *C, int size)
{
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            __m256 sum_vec = _mm256_setzero_ps();
            for(int k = 0; k < size; k+=8) {
                __m256 a_vec = _mm256_loadu_ps(A + i*size + k);
                __m256 b_vec = _mm256_loadu_ps(B + k*size + j);
                sum_vec = _mm256_add_ps(sum_vec, _mm256_mul_ps(a_vec, b_vec));
            }
            float sum = 0.0f;
            sum += sum_vec[0] + sum_vec[1] + sum_vec[2] + sum_vec[3]
                    + sum_vec[4] + sum_vec[5] + sum_vec[6] + sum_vec[7];
            C[i*size + j] = sum;
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, size;
    float *A, *B, *C, *buffer;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate matrices
    A = (float*) malloc(N*N*sizeof(float));
    B = (float*) malloc(N*N*sizeof(float));
    C = (float*) malloc(N*N*sizeof(float));

    // Initialize matrices
    for(int i = 0; i < N*N; i++) {
        A[i] = (float) rand()/RAND_MAX;
        B[i] = (float) rand()/RAND_MAX;
        C[i] = 0.0f;
    }

    // Scatter matrix A to all processes
    int rows_per_process = N/size;
    buffer = (float*) malloc(rows_per_process*N*sizeof(float));
    MPI_Scatter(A, rows_per_process*N, MPI_FLOAT, buffer, rows_per_process*N, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Compute matrix multiplication
    matrix_multiply(buffer, B, C, rows_per_process);

    // Gather results from all processes
    MPI_Gather(C, rows_per_process*N, MPI_FLOAT, C, rows_per_process*N, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Free memory
    free(A);
    free(B);
    free(C);
    free(buffer);

    MPI_Finalize();

    return 0;
}

//Part 2: Compare Performance of Vectorized and Non-Vectorized Matrix Multiplication

#define MAX_SIZE 2000

void matrix_multiply_naive(float *A, float *B, float *C, int size)
{
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                C[i*size + j] += A[i*size + k]*B[k*size + j];
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int rank, size, n_runs, n_sizes;
    float *A, *B, *C, *buffer;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0) {
        // Read command line arguments
        if(argc != 3) {
            printf("Usage: %s n_runs n_sizes\n", argv[0]);
            return 1;
        }
        n_runs = atoi(argv[1]);
        n_sizes = atoi(argv[2]);

        // Allocate matrices
        A = (float*) malloc(MAX_SIZE*MAX_SIZE*sizeof(float));
        B = (float*) malloc(MAX_SIZE*MAX_SIZE*sizeof(float));
        C = (float*) malloc(MAX_SIZE*MAX_SIZE*sizeof(float));

        // Initialize random seed
        srand(time(NULL));

        // Generate matrices of different sizes
        int sizes[n_sizes];
        for(int i = 0; i < n_sizes; i++) {
            sizes[i] = (i+1)*100;
        }

        // Print header
        printf("Size,Naive time,Vectorized time,Speedup\n");

        // Test matrix multiplication for each size
        for(int i = 0; i < n_sizes; i++) {
            int size = sizes[i];
            double naive_time = 0.0, vectorized_time = 0.0;

            // Perform test runs
            for(int j = 0; j < n_runs; j++) {
                // Initialize matrices
                for(int k = 0; k < size*size; k++) {
                    A[k] = (float) rand()/RAND_MAX;
                    B[k] = (float) rand()/RAND_MAX;
                }
                for(int k = 0; k < size*size; k++) {
                    C[k] = 0.0f;
                }

                // Perform naive multiplication
                double start_time = MPI_Wtime();
                matrix_multiply_naive(A, B, C, size);
                double end_time = MPI_Wtime();
                naive_time += end_time - start_time;

                // Perform vectorized multiplication
                for(int k = 0; k < size*size; k++) {
                    C[k] = 0.0f;
                }
                start_time = MPI_Wtime();
                matrix_multiply(A, B, C, size);
                end_time = MPI_Wtime();
                vectorized_time += end_time - start_time;
            }

            // Print performance results
            naive_time /= n_runs;
            vectorized_time /= n_runs;
            printf("%d,%.4f,%.4f,%.2f\n", size, naive_time, vectorized_time, naive_time/vectorized_time);
        }

        // Free matrices
        free(A);
        free(B);
        free(C);
    }

    MPI_Finalize();

    return 0;
}

