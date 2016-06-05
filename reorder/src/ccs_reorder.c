#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "papi.h"

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
    int *ccsRow, *ccsCol;
    double *value;
    int i, j;
    int cnt;
    double *vec, *b, *b_star;
    srand(time(NULL));
	float real_time, proc_time, mflops;
	long long flpins;
    
   
	// Read in the reordered(based on column) matrix
	FILE* pf;
    pf = fopen("matrix.reorder.ccs", "r");
    if (pf == NULL) {
        fprintf(stderr,  "Can't open input file!\n");
        return EXIT_FAILURE;
    }
    
	// Read in the matrix in AIJ Format
    fscanf(pf, "%d %d %d", &nRow, &nCol, &nnZero);
 
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

    // Allocate memory for CCS format
    ccsRow = (int*)malloc(sizeof(int)*nnZero);
    ccsCol = (int*)malloc(sizeof(int)*(nCol+1));
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in CCS!\n", nnZero, nCol+1, nnZero);
    memcpy(ccsRow, aijRow, nnZero*sizeof(int)); 
    cnt = 0;
    for (i = 0; i < nCol; i++) {
        ccsCol[i] = cnt;
        for (j = 0; j < nnZero; j++) {
            if (aijCol[j] == i)
                cnt++;
        }
    }
    ccsCol[i] = cnt;

    // Perform matrix-vector product
    // xA = b, A(nRow*nCol), x(1*nRow), b(1*nCol)
    vec = (double*)malloc(sizeof(double)*nRow);
    b = (double*)malloc(sizeof(double)*nCol);
    for(i = 0; i < nRow; ++i)
        vec[i] = rand()%100;
    for(i = 0; i < nCol; ++i)
        b[i] = 0.0;
    
    // Report Mflops/s rate using PAPI
    // Setup PAPI library and begin collecting array from the cnters
    PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    for (i = 0; i < nCol; i++) {
        b[i] = 0;
        for (j = ccsCol[i]; j < ccsCol[i + 1]; j++) {
            b[i] += vec[ccsRow[j]] * value[j];
        }
    }
    PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    printf("In reordered CCS format:\nReal_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n", real_time,  proc_time,  flpins,  mflops);
    PAPI_shutdown();

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

    b_star = (double*)malloc(sizeof(double)*nCol);
    for(i = 0; i < nCol; ++i)
        b_star[i] = 0.0;
    
    for(i = 0; i < nCol; ++i)
        for (j = 0; j < nRow; ++j)
            b_star[i] += matrixA[j][i] * vec[j];
    
    if(norm(b, b_star, nCol) < 1e-6)
    {
        printf("The program of vec-mat product is right!\n");
    }
    else
    {
        printf("The error is too large and valued as %lf!\n", norm(b, b_star, nCol));
    }
    return EXIT_SUCCESS;
}
