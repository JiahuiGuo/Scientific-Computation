/**********************************************
 * Author: Jiahui Guo
 * Email: guo.jiahui07@gmail.com
 * Department: EECS
 * CS594 Homework 1
 * Function: The matrx-vector multiplication
***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "../include/c_timer.h"
#include "../include/cblas.h"

//#define DEBUG
#define MAX 100
#define MIN 0

void vecMul(double *vec_y, double *vec_x, double **matrixA, size_t size)
{
	size_t i, j;
    for(i = 0; i < size; ++i)
        for(j = 0; j < size; ++j)
            vec_y[i] += matrixA[i][j] * vec_x[j];
}

bool verify(double *vec_y, double *vec_x, double **matrixA, size_t size)
{
    bool status = false;
    double *aVec_y;
    double *error;
    aVec_y = (double*)malloc(size*sizeof(double));
    error = (double*)malloc(size*sizeof(double));
    //void    cblas_dgemv (const enum CBLAS_ORDER order,  const enum CBLAS_TRANSPOSE TransA,  const int M,  const int N,  const double alpha,  const double *A,  const int lda,  const double *X,  const int incX,  const double beta,  double *Y,  const int incY)
    //        Multiplies a matrix and a vector. )
    double *mA;
    mA = (double*)calloc((size*size), sizeof(double));
    //Convert 2D matrixA to 1D mA
	size_t i, j, k;
    for(i = 0; i < size; ++i)
        for(j = 0; j < size; ++j)
            mA[i*size + j] = matrixA[i][j];
    //Calculate the result using a standard lib
    cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1, mA, size, vec_x, 1, 1, aVec_y, 1);
    //Obtain the error terms
    for(k = 0; k < size; ++k)
        error[k] = fabs(aVec_y[k] - vec_y[k]); 
    //Check whether the result is accurate
    status = (cblas_dnrm2(size, error, 1) < 1e-16);
    //Free the memory
    free(aVec_y);
    free(error);
    aVec_y = NULL;
    error = NULL;

#ifdef DEBUG
    printf("The verification result is %d\n", status);
#endif
    
    return status;
}

int main(int argc, char* argv[])
{
    size_t size;
    double btime, etime;
    double *vec_x;
    double *vec_y;
    double **matrixA;

    srand(time(NULL));
    if(argc > 1)
    {
        size = atoi(argv[1]);
#ifndef DEBUG
        printf("%zd\t", size);
#endif
    }
    else
    {
        printf("Not enough parameters!\n");
        return EXIT_FAILURE;
    }
    
    vec_x = (double*)malloc(size*sizeof(double));
    vec_y = (double*)malloc(size*sizeof(double));
    matrixA = (double**)malloc(size*sizeof(double));
    size_t i;
	for(i = 0; i < size; ++i)
        matrixA[i] = (double*)malloc(size*sizeof(double));

    if(vec_x == NULL || vec_y == NULL || matrixA == NULL)
    {
        printf("The memory allocation failed! Exit Now!\n");
        return EXIT_FAILURE;
    }
    else
    {
        //initialize the array
        size_t j;
		for(j = 0; j < size; ++j)
        {
            vec_x[j] = rand()%(MAX-MIN+1);
            vec_y[j] = 0.0;
        }
        //initialize the matrix
        for(i = 0; i < size; ++i)
            for(j = 0; j < size; ++j)
                matrixA[i][j] = rand()%(MAX-MIN+1);

    }

    //Calculate the execution time of vector and matrix multiplication
    btime = get_cur_time();
    vecMul(vec_y, vec_x, matrixA, size);
    etime = get_cur_time();

    //verify the result
    if(verify(vec_y, vec_x, matrixA, size))
    {
#ifdef DEBUG
        printf("The result is accurate!\n");
#endif
    }
    else
    {
#ifdef DEBUG
        printf("The result is not acccurate! Exit Now!\n");
#endif
        return EXIT_FAILURE;
    }

    //Free the allocated memory
    free(vec_x);
    free(vec_y);
    vec_x = NULL;
    vec_y = NULL;
    for(i = 0; i < size; ++i)
    {
        free(matrixA[i]);
        matrixA[i] = NULL;
    }   
    free(matrixA);
    matrixA = NULL;

#ifdef DEBUG
    printf("Elapsed time is %lf seconds\n", etime - btime);
    printf("The FLOPS is %lf GOps/sec\n", size * size * 2 / (etime - btime) / 1e+9);
#else
    printf("%.16f\t", etime - btime);
    printf("%.16f\n", size * size * 2 / (etime - btime) / 1e+9);
#endif

    return EXIT_SUCCESS;
}
