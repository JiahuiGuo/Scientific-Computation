/**********************************************
 * Author: Jiahui Guo
 * Email: guo.jiahui07@gmail.com
 * Department: EECS
 * CS594 Homework 1
 * Function: The 2-norm of a vector
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

double twoNorm(const double *data, const size_t size)
{
    double sum = 0.0;
    size_t i;
    for(i = 0; i < size; ++i)
    {
        sum += data[i] * data[i];
    }
#ifdef DEBUG
    printf("The result is %.4f\n", sqrt(sum));
#endif
    return sqrt(sum);
}

bool verify(const double result, const double *data, const size_t size)
{
    bool status = false;
    // double cblas_dnrm2(const int N,  const double *X,  const int incX);
    double stdResult = cblas_dnrm2(size, data, 1);
    status = (abs(stdResult - result) < 1e-16);

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
    double result;
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
    if(vec_x == NULL)
    {
        printf("The memory allocation failed! Exit Now!\n");
        return EXIT_FAILURE;
    }
    else
    {
        size_t i;
        for(i = 0; i < size; ++i)
        {
            vec_x[i] = rand()%(MAX-MIN+1);
#ifdef DEBUG
            printf("The random data generated is %.4f\n", vec_x[i]);
#endif
        }
    }

    //Calculate the execution time of the 2 norm of a vector
    btime = get_cur_time();
    result = twoNorm(vec_x, size);
    etime = get_cur_time();

    //verify the result
    if(verify(result, vec_x, size))
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
    vec_x = NULL;

#ifdef DEBUG
    printf("Elapsed time is %lf seconds\n", etime - btime);
    printf("The FLOPS is %lf GOps/sec\n", size * 2 / (etime - btime));
#else
    printf("%.16f\t", etime - btime);
    printf("%.16f\n", size * 2 / (etime - btime) / 1e+9);
#endif

    return EXIT_SUCCESS;
}
