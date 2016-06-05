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
    int *crsRow, *crsCol;
    int *bcrsRow, *bcrsCol;
    double *value;
    int i, j, k;
    int cnt;
    double *vec, *b_star;
    double *bcrsData, *bcrsResult;
	int sub, index, flag;
	double *blkinit;
    srand(time(NULL));
    float real_time, proc_time, mflops;
	long long flpins;
    
    // Open the file
    FILE* pf;
    pf = fopen("matrix.reorder.crs", "r");
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

    // Allocate memory for CRS format
    crsRow = (int*)malloc(sizeof(int)*(nRow+1));
    crsCol = (int*)malloc(sizeof(int)*nnZero);
    printf("The matrix is has %d(Row), %d(Col), %d(non-zero) in CRS!\n", nRow+1, nnZero, nnZero);
    memcpy(crsCol, aijCol, nnZero*sizeof(int)); 
    cnt = 0;
    for (i = 0; i < nRow; ++i) {
        crsRow[i] = cnt;
        for (j = 0; j < nnZero; ++j) {
            if (aijRow[j] == i)
                cnt++;
        }
    }
    crsRow[i] = cnt;

    bcrsRow = (int*)malloc(sizeof(int)*(nRow/3 + 1));
    bcrsCol = (int*)malloc(sizeof(int)*(nRow/3)*(nRow/3));
    bcrsData = (double*)malloc(sizeof(double)*nnZero*9);
    bcrsResult = (double*)malloc(sizeof(double)*nRow);

	for (i = 0; i < nnZero*9; ++i) {
		bcrsData[i] = 0;
	}
	
    vec = (double*)malloc(sizeof(double)*nRow);
    for(i = 0; i < nRow; ++i)
        vec[i] = rand()%100;

    // BCRS
	cnt = 0;
	for (i = 0; i < nRow/3; ++i) { //ranging over bcrsRow
		bcrsRow[i] = cnt;
		for (sub = 3 * i; sub < 3 * i + 3; ++sub) { // sub ranging over crsRow
			for (j = crsRow[sub]; j < crsRow[sub + 1]; ++j) { //j ranging over crsCol
				blkinit = bcrsData + 9 * bcrsRow[i]; // blkinit -- the starting point of the block correspoinding to i
				flag = 0;
				for (index = bcrsRow[i]; index < cnt; ++index) {
					if (crsCol[j]/3 == bcrsCol[index]) {
						flag = 1;
						break;
					}
					blkinit = bcrsData + 9 * (index + 1);
				}
				*(blkinit + (crsCol[j]%3) + (sub%3)*3) = value[j];
				if (flag == 0) {
					bcrsCol[cnt] = crsCol[j]/3;
					cnt++;
				}
			}
		}
	}
	bcrsRow[i] = cnt;
	
	PAPI_flops(&real_time, &proc_time, &flpins, &mflops);
	for (i = 0; i < nRow/3; ++i) {
		*(bcrsResult + 3*i) = 0;
		*(bcrsResult + 3*i+1) = 0;
		*(bcrsResult + 3*i+2) = 0;
		for (j = bcrsRow[i]; j < bcrsRow[i + 1]; ++j) {
			for (k = 0; k < 9; ++k) {
				*(bcrsResult + 3*i + k/3) += bcrsData[9*j + k] * vec[bcrsCol[j]*3 + k%3];
			}
		} 
	}
	PAPI_flops(&real_time, &proc_time, &flpins, &mflops);
	printf("In BCRS format:\nReal_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\nBRCS blocks:\t%d\n", real_time, proc_time, flpins, mflops, cnt);
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

    b_star = (double*)malloc(sizeof(double)*nRow);
    for(i = 0; i < nRow; ++i)
        b_star[i] = 0.0;
    
    for(i = 0; i < nRow; ++i)
        for (j = 0; j < nCol; ++j)
            b_star[i] += matrixA[i][j] * vec[j];
    
    if(norm(bcrsResult, b_star, nRow) < 1e-6)
    {
        printf("The program of vec-mat product is right!\n");
    }
    else
    {
        printf("The error is too large and valued as %lf!\n", norm(bcrsResult, b_star, nRow));
    }
    return EXIT_SUCCESS;
}
