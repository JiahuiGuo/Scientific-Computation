#ifndef _PDGEMM_H_
#define _PDGEMM_H_

#define NUM_THREADS 8
#define MAX 10
#define MIN 0

typedef struct {
	int row;
	int col;
	int b_row;			// row of blocks per process
	int b_col;			// column of blocks per process
	int ld_row;			// row of total elements per process
	int ld_col; 		// column of total elements per process
	double *mat;
} proc_mat;


typedef struct {
	int id;
	int thread_row;
	int thread_col;
	double* thread_A;
	double* thread_B;
	double* thread_C;
} thread_data;

void InitMat(double *, int, int);
void PrintMat(double *, int);
void findBlock(proc_mat*, int, int, proc_mat*, int, int);
void cyclicData(proc_mat*, proc_mat*);
void colToBuf(proc_mat*, double*, int, int);
void rowToBuf(proc_mat*, double*, int, int);
void collect(int, proc_mat*, proc_mat*, int);
void* mmm(void*);
void pdgemm(double*, double*, double*, int);

#endif
