#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
//#include "papi.h"
//

// This function is used to count the number of the identical menber of one array 
int getCompressed(int N, int* array)
{
    int size = 1;
    int i, tmp;
    tmp = array[0];
    for(i=1; i < N; ++i){
        if(tmp != array[i]){
            size += 1;
            tmp = array[i];
        }
    }
    return size + 1;
}

// This fucntion is used to copy the identical member of one array to another
// array1 <- array2
void getValue(int* array1, int* array2, int N)
{
    int i, tmp;
    int j = 0;
    tmp = array2[0];
    array1[0] = array2[0];
    for(i = 0; i < N; ++i)
    {
        if(tmp != array2[i])
        {
            tmp = array2[i];
            j++;
            array1[j] = i + 1;
        }
    }
    array1[j+1] = N + 1;
}

double norm(const double *data1, const double *data2, const size_t size)
{
    double* error;
    error = (double*)malloc(sizeof(double)*size);
    double sum = 0.0;
    size_t i;
    for(i = 0; i < size; ++i)
    {
        error[i] = fabs(data1[i] - data2[i]);
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
    int nCrsRow = 0;
    double *value;
    int i, j, k;
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
    }
    fclose(pf);

    // It is impossible to get the length of a dynamic allocated array using sizeof!!
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in AIJ!\n", nnZero, nnZero, nnZero);

    // Allocate memory for CRS format
    nCrsRow = getCompressed(nnZero, aijRow);
    crsRow = (int*)malloc(sizeof(int)*nCrsRow);
    crsCol = (int*)malloc(sizeof(int)*nnZero);
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in CRS!\n", nCrsRow, nnZero, nnZero);
    getValue(crsRow, aijRow, nnZero);    
    memcpy(crsCol, aijCol, nnZero*sizeof(int));
    
    printf("\ncrsRow is: ");
    for(i = 0; i < nCrsRow; ++i)
        printf("%d\t", crsRow[i]);
    printf("\ncrsCol is: ");
    for(i = 0; i < nnZero; ++i)
        printf("%d\t", crsCol[i]);
  

    // Write the CRS format to a file
    pf = fopen("matrix.output.crs", "w");
    if (pf == NULL)
    {
        fprintf(stderr, "Can't open output file!\n");
        return EXIT_FAILURE;
    }

    k = 0;
    fprintf(pf, "%d %d %d\n", nCrsRow, nnZero, nnZero);
    for(i = 0; i < nCrsRow; ++i)
    {
        fprintf(pf, "%d\t", crsRow[i]);
        for(j = crsRow[i]; j < crsRow[i+1]; ++j)
        {
            fprintf(pf, "\t\t%d\t%lf\n", crsCol[k], value[k]);
            k++;
        }
    }
    

    // Perform matrix-vector product
    // Ax = b, A(nRow*nCol), x(nCol*1), b(nRow*1)
    vec = (double*)malloc(sizeof(double)*nCol);
    b = (double*)malloc(sizeof(double)*nRow);
    for(i = 0; i < nCol; ++i)
        vec[i] = rand()%100;
    for(i = 0; i < nRow; ++i)
        b[i] = 0.0;
    
    // Report Mflops/s rate using PAPI
    /* Setup PAPI library and begin collecting data from the counters */
    //PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    
    k = 0;
    for(i = 0; i < nCrsRow; ++i)
        for(j = crsRow[i]; j < crsRow[i+1]; ++j)
        {
            b[i] +=  value[k] * vec[crsCol[j-1]-1];
            k++;
        }
   
    //PAPI_flops(&real_time,  &proc_time,  &flpins,  &mflops);
    //printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n", real_time,  proc_time,  flpins,  mflops);
    //PAPI_shutdown();
     
    printf("vec:\n");
    for (i = 0; i < nCol; ++i)
        printf("%lf\t", vec[i]);
    
    printf("\nb:\n");
    for (i = 0; i < nRow; ++i)
        printf("%lf\t", b[i]);

    // Verify the correctness of the code
    // Construct the original matrix
    matrixA = (double**)malloc(sizeof(double)*nRow);
    for(i = 0; i < nRow; ++i)
        matrixA[i] = (double*)malloc(sizeof(double)*nCol);
    for(i = 0; i < nRow; ++i)
        for(j = 0; j < nCol; ++j)
            matrixA[i][j] = 0.0;
    for(i = 0; i < nnZero; ++i)
        matrixA[aijRow[i]-1][aijCol[i]-1] = value[i];
    
    for(i = 0; i < nRow; ++i)
    {
        for(j = 0; j < nCol; ++j)
            printf("%lf\t", matrixA[i][j]);
        printf("\n");
    }
    b_star = (double*)malloc(sizeof(double)*nCol);
    for(i = 0; i < nRow; ++i)
        b_star[i] = 0.0;
    for(i = 0; i < nRow; ++i)
        for (j = 0; j < nCol; ++j)
            b_star[i] += matrixA[i][j] * vec[j];
    
    printf("\nb_star:\n");
    for (i = 0; i < nRow; ++i)
        printf("%lf\t", b_star[i]);
   
   
    if(norm(b, b_star, nCol) < 1e-16)
    {
        printf("The program of vec-mat product is right!\n");
    }
    else
    {
        printf("The error is too large and valued as \n!%lf", norm(b, b_star, nCol));
    }
    return EXIT_SUCCESS;
}
