#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
//#include "papi.h"
//

double norm(const double *array1, const double *array2, const size_t size)
{
    double* error;
    error = (double*)malloc(sizeof(double)*size);
    double sum = 0.0;
    size_t i;
    for(i = 0; i < size; ++i)
    {
        error[i] = fabs(array1[i] - array2[i]);
        sum += error[i] * error[i];
    }
    free(error);
    return sqrt(sum);
}

int main()
{
    int nRow, nCol, nnZero;
    double** matrixA;
    int *aijRow, *aijCol;
    int *crsRow, *crsCol;
    double *value;
    int i, j;
    int cnt;
    double *vec, *b, *b_star;
    srand(time(NULL));
    // float real_time, proc_time, mflops;
    //long long flpins;
    
    // Open the file
    FILE* pf;
    //pf = fopen("matrix.output", "r");
    pf = fopen("matrix.test", "r");
    if (pf == NULL) {
        fprintf(stderr,  "Can't open input file!\n");
        return EXIT_FAILURE;
    }
    // Read in the matrix in AIJ Format
    fscanf(pf, "%d %d %d", &nRow, &nCol, &nnZero);
    printf("The matrix is in %d(Row) * %d(Col) size with %d non-zero elements!\n", nRow, nCol, nnZero);
    
    // Allocate memory for AIJ format
    aijRow = (int*)malloc(sizeof(int)*nnZero);
    aijCol = (int*)malloc(sizeof(int)*nnZero);
    value = (double*)malloc(sizeof(double)*nnZero);
    
    // Read in the AIJ format matrix
    for(i = 0; i < nnZero; ++i){
        // Note: The index stored is not C style!
        fscanf(pf, "%d %d %lf", &aijRow[i], &aijCol[i], &value[i]);   
        // Change to C style
        aijRow[i] -= 1;
        aijCol[i] -= 1;
    }
    fclose(pf);

    // It is impossible to get the length of a dynamic allocated array using sizeof!!
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in AIJ!\n", nnZero, nnZero, nnZero);

    // Sort the AIJ format with the row indices in ascending order

    // Allocate memory for CRS format
    crsRow = (int*)malloc(sizeof(int)*(nRow+1));
    crsCol = (int*)malloc(sizeof(int)*nnZero);
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in CRS!\n", nRow + 1, nnZero, nnZero);
    memcpy(crsCol, aijCol, nnZero*sizeof(int)); 
    cnt = 0;
    for (i = 0; i < nRow; i++) {
        crsRow[i] = cnt;
        for (j = 0; j < nnZero; j++) {
            if (aijRow[j] == i)
                cnt++;
        }
    }
    crsRow[i] = cnt;

    // Perform matrix-vector product
    // Ax = b, A(nRow*nCol), x(nCol*1), b(nRow*1)
    vec = (double*)malloc(sizeof(double)*nCol);
    b = (double*)malloc(sizeof(double)*nRow);
    for(i = 0; i < nCol; ++i)
        vec[i] = rand()%100;
    for(i = 0; i < nRow; ++i)
        b[i] = 0.0;
    
    // Report Mflops/s rate using PAPI
    /* Setup PAPI library and begin collecting array from the cnters */
    //PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    for (i = 0; i < nRow; i++) {
        b[i] = 0;
        for (j = crsRow[i]; j < crsRow[i + 1]; j++) {
            b[i] += value[j]*vec[crsCol[j]];
        }
    }
   
    //PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    //printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n", real_time,  proc_time,  flpins,  mflops);
    //PAPI_shutdown();
     
    // Verify the correctness of the code
    // Construct the original matrix
    matrixA = (double**)malloc(sizeof(double)*nRow);
    for(i = 0; i < nRow; ++i)
        matrixA[i] = (double*)malloc(sizeof(double)*nCol);
    for(i = 0; i < nRow; ++i)
        for(j = 0; j < nCol; ++j)
            matrixA[i][j] = 0.0;
    
    for(i = 0; i < nnZero; ++i)
        matrixA[aijRow[i]][aijCol[i]] = value[i];
    /*
    for(i = 0; i < nRow; ++i)
    {
        for(j = 0; j < nCol; ++j)
            printf("%lf\t", matrixA[i][j]);
        printf("\n");
    }*/

    b_star = (double*)malloc(sizeof(double)*nRow);
    for(i = 0; i < nRow; ++i)
        b_star[i] = 0.0;
    for(i = 0; i < nRow; ++i)
        for (j = 0; j < nCol; ++j)
            b_star[i] += matrixA[i][j] * vec[j];
    
    if(norm(b, b_star, nCol) < 1e-16)
    {
        printf("The program of vec-mat product is right!\n");
    }
    else
    {
        printf("The error is too large and valued as %lf\n!", norm(b, b_star, nCol));
    }
    return EXIT_SUCCESS;
}
