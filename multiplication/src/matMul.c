/**********************************************
 * Function: The matrix-matrix multiplication
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

void matMul(double **matrixA, double **matrixB, double **matrixC, size_t size)
{	
	size_t i, j, k;
    for(i = 0; i < size; ++i)
        for(j = 0; j < size; ++j)
            for(k = 0; k < size; ++k)
                matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
}

bool verify(double **matrixA, double **matrixB, double **matrixC, size_t size)
{
    bool status = false;
    double *aMatrixC;
    double *error;
    
    aMatrixC = (double*)malloc(size*size*sizeof(double));
    error = (double*)malloc(size*size*sizeof(double));
    
    double *mA;
    mA = (double*)calloc((size*size), sizeof(double));
    double *mB;
    mB = (double*)calloc((size*size), sizeof(double));
    //Convert 2D matrixA to 1D mA
    size_t i, j;
	for(i = 0; i < size; ++i)
    {   
        for(j = 0; j < size; ++j)
        {
            mA[i*size + j] = matrixA[i][j];
            mB[i*size + j] = matrixB[i][j];
        }
    }
    //Calculate the result using a standard lib
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, mA, size, mB, size, 1, aMatrixC, size);
    //Obtain the error terms
    for(i = 0; i < size; ++i)
        for(j = 0; j < size; ++j)
            error[i*size + j] = fabs(aMatrixC[i*size + j] - matrixC[i][j]);
    //Check whether the result is accurate
    status = (cblas_dnrm2(size * size, error, 1) < 1e-16);
    //Free the memory
    free(aMatrixC);
    aMatrixC = NULL;
    free(error);
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
    double **matrixA;
    double **matrixB;
    double **matrixC;

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
    
	size_t i;
    matrixA = (double**)malloc(size*sizeof(double));
    for(i = 0; i < size; ++i)
        matrixA[i] = (double*)malloc(size*sizeof(double));
 
    matrixB = (double**)malloc(size*sizeof(double));
    for(i = 0; i < size; ++i)
        matrixB[i] = (double*)malloc(size*sizeof(double));

   matrixC = (double**)malloc(size*sizeof(double));
    for(i = 0; i < size; ++i)
        matrixC[i] = (double*)malloc(size*sizeof(double));

    if(matrixA == NULL || matrixB == NULL || matrixC == NULL)
    {
        printf("The memory allocation failed! Exit Now!\n");
        return EXIT_FAILURE;
    }
    else
    {
		size_t j;
        //initialize the matrix
        for(i = 0; i < size; ++i)
        {
            for(j = 0; j < size; ++j)
            {
                matrixA[i][j] = rand()%(MAX-MIN+1);
                matrixB[i][j] = rand()%(MAX-MIN+1);
                matrixC[i][j] = 0.0;
            }
        }
    }

    //Calculate the execution time of the 2 norm of a vector
    btime = get_cur_time();
    matMul(matrixA, matrixB, matrixC, size);
    etime = get_cur_time();

    //verify the result
    if(verify(matrixA, matrixB, matrixC, size))
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
    for(i = 0; i < size; ++i)
    {
        free(matrixA[i]);
        free(matrixB[i]);
        free(matrixC[i]);
        matrixA[i] = NULL;
        matrixB[i] = NULL;
        matrixC[i] = NULL;
    }   
    free(matrixA);
    free(matrixB);
    free(matrixC);
    matrixA = NULL;
    matrixB = NULL;
    matrixC = NULL;

#ifdef DEBUG
    printf("Elapsed time is %lf seconds\n", etime - btime);
    printf("The FLOPS is %lf GOps/sec\n", size * size * size * 2 / (etime - btime) / 1e+9);
#else
    printf("%.16f\t", etime - btime);
    printf("%.16f\n", size * size * size * 2 / (etime - btime) / 1e+9) ;
#endif

    return EXIT_SUCCESS;
}
