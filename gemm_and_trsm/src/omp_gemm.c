/**********************************************
 * Author: Jiahui Guo
 * Email: guo.jiahui07@gmail.com
 * Department: EECS
 * CS594 Homework 2
 * Function: Implementation of general matrix-matrix
 *           multiplication in OpenMP
***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "../include/c_timer.h"

//#define DEBUG
#define MAX 100
#define MIN 0
#define CHUNKSIZE 10
#define THREADS 2

#define MM 400
#define NN 400
#define KK 400
    
//gemm: C <- alpha*A*B + beta*C
#define ALPHA -1
#define BETA 1

int main(int argc, char* argv[])
{
    int tid, chunk, nthreads;
	double g_btime, g_etime;
    size_t M, N, K;
    double **matrixA;
    double **matrixB;
    double **matrixC;

    srand(time(NULL));

    if(argc > 3)
    {
        M = atoi(argv[1]);
        N = atoi(argv[2]);
        K = atoi(argv[3]);
#ifdef DEBUG
        printf("The demensions for hte matrices are M=%zd, N=%zd, K=%zd.\n", M, N, K);
#endif
    }
    else
    {
        M = MM;
        N = NN;
        K = KK;
#ifdef DEBUG
        printf("Use the default parameters for GEMM!\n");
#endif
    }
    
    chunk = CHUNKSIZE;

	size_t m, n, k;
    matrixA = (double**)malloc(M*sizeof(double));
    for(m = 0; m < M; ++m)
        matrixA[m] = (double*)malloc(N*sizeof(double));
 
    matrixB = (double**)malloc(N*sizeof(double));
    for(n = 0; n < N; ++n)
        matrixB[n] = (double*)malloc(K*sizeof(double));

    matrixC = (double**)malloc(M*sizeof(double));
    for(m = 0; m < M; ++m)
        matrixC[m] = (double*)malloc(K*sizeof(double));

    if(matrixA == NULL || matrixB == NULL || matrixC == NULL)
    {
        printf("The memory allocation failed! Exit Now!\n");
        return EXIT_FAILURE;
    }
    
    omp_set_num_threads(THREADS);
    g_btime = get_cur_time();

    #pragma omp parallel shared(matrixA, matrixB, matrixC, nthreads, chunk) private(tid, m, n, k)
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
        
        #pragma omp for schedule(static, chunk)
        for(m = 0; m < M; ++m)
            for(n = 0; n < N; ++n)
                matrixA[m][n] = rand()%(MAX-MIN+1);
        #pragma omp for schedule(static, chunk)
        for(n = 0; n < N; ++n)
            for(k = 0; k < K; ++k)
                matrixB[n][k] = rand()%(MAX-MIN+1);
        #pragma omp for schedule(static, chunk)
        for(m = 0; m < M; ++m)
            for(k = 0; k < K; ++k)
                matrixC[m][k] = 0.0;
        #pragma omp for schedule(static, chunk)
        for(m = 0; m < M; ++m)
            for(n = 0; n < N; ++n)
                for(k = 0; k < K; ++k)
                    matrixC[m][k] = ALPHA * matrixA[m][n] * matrixB[n][k] + BETA * matrixC[m][k];
    } //End of parallel section
    g_etime = get_cur_time();

#ifdef DEBUG
    printf("Elapsed time for the whole program is %lf seconds.\n", g_etime - g_btime);
	printf("Matrix A is:\n");
    for(m = 0; m < M; ++m)
    {
		for(n = 0; n < N; ++n)
			printf("%f\t", matrixA[m][n]);
        printf("\n");
    }

	printf("Matrix B is:\n");
    for(n = 0; n < N; ++n)
    {
        for(k = 0; k < K; ++k)
            printf("%f\t", matrixB[n][k]);
        printf("\n");
    }

	printf("Matrix C is:\n");
    for(m = 0; m < M; ++m)
    {
        for(k = 0; k < K; ++k)
            printf("%f\t", matrixC[m][k]);
        printf("\n");
    }
#else
        printf("%zd\t%.16f\n", M, g_etime - g_btime);
#endif

    //Free the allocated memory
    for(m = 0; m < M; ++m)
    {    
        free(matrixA[m]);
        matrixA[m] = NULL;
    }
    for(n = 0; n < N; ++n)
    {   
        free(matrixB[n]);
        matrixB[n] = NULL;
    }
    for(m = 0; m < M; ++m)
    {
        free(matrixC[m]);
        matrixC[m] = NULL;
    }
    free(matrixA);
    free(matrixB);
    free(matrixC);
    matrixA = NULL;
    matrixB = NULL;
    matrixC = NULL;

    return EXIT_SUCCESS;
}
