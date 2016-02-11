/****************************************************
 * Author: Jiahui Guo
 * Email: guo.jiahui07@gmail.com
 * Department: EECS
 * CS594 Homework 2
 * Function: Implementation of Triangular Solve with 
 *           Matrix in OpenMP
 * Description: The TRSM routines solves AX=B system
 *  where A is either upper or lower triangular. The
 *  solution matrix X is returned in the space occupied
 *  by B
****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include "../include/c_timer.h"

//#define DEBUG
#define MAX 100
#define MIN 0
#define CHUNKSIZE 10
#define THREADS 8

#define MM 400
#define NN 400

//trsm: AX = B

int main(int argc, char* argv[])
{
    int tid, chunk, nthreads;
    double g_btime, g_etime;
    size_t M, N;
    size_t k;
    double **matrixA;
    double **matrixB;

    srand(time(NULL));

    if(argc > 2)
    {
        M = atoi(argv[1]);
        N = atoi(argv[2]);
#ifdef DEBUG
        printf("The demensions for hte matrices are M=%zd, N=%zd.\n", M, N);
#endif
    }
    else
    {
        M = MM;
        N = NN;
#ifdef DEBUG
        printf("Use the default parameters for GEMM!\n");
#endif
    }
    
    chunk = CHUNKSIZE;

	size_t i, j;
    matrixA = (double**)malloc(M*sizeof(double));
    for(i = 0; i < M; ++i)
        matrixA[i] = (double*)malloc(M*sizeof(double));
 
    matrixB = (double**)malloc(M*sizeof(double));
    for(j = 0; j < M; ++j)
        matrixB[j] = (double*)malloc(N*sizeof(double));

    if(matrixA == NULL || matrixB == NULL)
    {
        printf("The memory allocation failed! Exit Now!\n");
        return EXIT_FAILURE;
    }
    

    omp_set_num_threads(THREADS);
    g_btime = get_cur_time();
    #pragma omp parallel shared(matrixA, matrixB, nthreads, chunk) private(tid, i, j, k)
    { 
        tid = omp_get_thread_num();
        if(tid == 0)    //master thread
        {
            nthreads = omp_get_num_threads();
#ifdef DEBUG
            printf("Solving GEMM using %d threads!\n", nthreads);
            printf("The number of processors available is %d.\n", omp_get_num_procs());
            printf("The number of threads been used is %d.\n", omp_get_num_threads());
#endif
        }
    
	// Initialization
    // Matrix A is a lower triangular
    #pragma omp for schedule(static, chunk)
	for(i = 0; i < M; ++i)
        for(j = 0; j < M; ++j)
            if(i >= j)
                matrixA[i][j] = rand()%(MAX-MIN+1);
            else
                matrixA[i][j] = 0;

    // Matrix B is full matrix
    #pragma omp for schedule(static, chunk)
	for(i = 0; i < M; ++i)
        for(j = 0; j < N; ++j)
            matrixB[i][j] = rand()%(MAX-MIN+1);
	
	//Solve the equation
	#pragma omp for schedule(static, chunk)	
    for(j = 0; j < N; j++)
    {
        for(k = 0; k < M; k++)
        {
            matrixB[k][j] = matrixB[k][j] / matrixA[k][k];
            for(i = k + 1; i < M; i++)
                matrixB[i][j] = matrixB[i][j] - matrixB[k][j] * matrixA[i][k];
        }
  	}
	
	}
    //End of parallel section
 
    g_etime = get_cur_time();
#ifdef DEBUG
    printf("The matrix A is:\n");
    for(i = 0; i < M; ++i)
	{
        for(j = 0; j < M; ++j)
            printf("%f\t", matrixA[i][j]);
        printf("\n");
    }
    printf("The matrix B is:\n");
    for(i = 0; i < M; ++i)
    {
        for(j = 0; j < N; ++j)
            printf("%f\t", matrixB[i][j]);
        printf("\n");
    }
#else
    printf("%zd\t%.16f\n", M, g_etime - g_btime);
#endif

    //Free the allocated memory
    for(i = 0; i < M; ++i)
    {    
        free(matrixA[i]);
        matrixA[i] = NULL;
    }
    for(j = 0; j < M; ++j)
    {   
        free(matrixB[j]);
        matrixB[j] = NULL;
    }
    free(matrixA);
    free(matrixB);
    matrixA = NULL;
    matrixB = NULL;

    return EXIT_SUCCESS;
}
