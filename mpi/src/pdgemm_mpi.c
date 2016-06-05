/*
	The function implements a MPI based MMM
	Assume that P = Q = sqrt(size), the matrices are all N-by-N, 
	and N is divisible by sqrt(size).
	Author: Jiahui Guo
	Email: jguo7@utk.edu

	Acknowlegement: Some of the ideas are from Dr.Bosilca, Chunyan, and Chongxiao.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <cblas.h>
#include <pthread.h>
#include "pdgemm_mpi.h"

int block;
int P, Q;
int rank, size;

int main(int argc, char** argv){

    int N;
	if(argc > 3)
	{
		P = atoi(argv[1]);
		Q = atoi(argv[2]);
		N = atoi(argv[3]);
	}
	else
	{		
		printf("Usage: ./pdgemm P Q N");
		exit(1);
	}
	block = N * N;
	double *A, *B, *C;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0){
		// Allocation all the matrices in processor 0
		A = (double*)malloc(N*N*sizeof(double));
		B = (double*)malloc(N*N*sizeof(double));
		C = (double*)malloc(N*N*sizeof(double));
		if(A == NULL || B == NULL || C == NULL)
        {
            printf("Memory is not allocated!\n");
            exit(1);
        }
        else{
            InitMat(A, N, 1);
		    InitMat(B, N, 1);
		    InitMat(C, N, 0);
        }
        printf("Matrix A is:\n");
        PrintMat(A, N);
        printf("Matrix B is:\n");
        PrintMat(B, N);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	pdgemm(A, B, C, N);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("Matrix C is:\n");
		PrintMat(C, N);
		free(A);
		free(B);
		free(C);
	}
	MPI_Finalize();
	return 0;
}

void InitMat(double *m, int n, int flag){
	// flag is indicated whether it is a matrix with all 0's
	srand((unsigned)time(NULL));
	int i;
	for(i = 0; i < n*n; ++i){
		if(flag)
			m[i] = rand()%(MAX-MIN+1);
		else
			m[i] = 0.0;
	}
}

void PrintMat(double *m, int n)
{
	size_t i, j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j)
			printf("%f\t", m[j + i * n]);
		printf("\n");
	}
	printf("\n");
}

void findBlock(proc_mat *proc_A, int local_r, int local_c, proc_mat *old_A, int index_r, int index_c)
{
	int r,c;
	int r_block = proc_A->ld_row/proc_A->b_row;
	int c_block = proc_A->ld_col/proc_A->b_col;
	int proc_i, proc_j, global_i, global_j;
	for(r = 0; r < r_block; ++r){
		for(c = 0; c < c_block; ++c){
			proc_i = local_r * r_block+r;
			proc_j = local_c * c_block+c;
			global_i = index_r * r_block+r;
			global_j = index_c * c_block+c;
			if((global_i <= old_A->ld_row-1) && (global_j <= old_A->ld_col-1)){
				proc_A->mat[proc_j*proc_A->ld_row + proc_i] = old_A->mat[global_j*old_A->ld_row + global_i];
				if(local_r == 0){
					if(proc_A->col < proc_j + 1){
						proc_A->col = proc_j + 1;
					}
				}
				if(local_c == 0){
					if(proc_A->row < proc_i + 1){
						proc_A->row = proc_i + 1;
					}
				}
			}
		}
	}
}

void cyclicData(proc_mat *proc_A, proc_mat *old_A){
	int i, id, index_r, index_c, local_r, local_c;
	for(id = size - 1; id >= 0; --id){
		for(i = 0; i < (proc_A->ld_row)*(proc_A->ld_col); ++i){
			proc_A->mat[i] = 0;
		}	
		local_c = 0;
		for(index_c = id%Q; index_c < old_A->b_col; index_c += Q){
			local_r = 0;
			for(index_r = id/Q; index_r < old_A->b_row; index_r += P){
				findBlock(proc_A, local_r, local_c, old_A, index_r, index_c);
				local_r++;
			}
			local_c++;
		}
		// Send proc_mat and data to each process
		// tag = 0
		if(id != 0){
			MPI_Send(proc_A, sizeof(proc_mat), MPI_BYTE, id, 0, MPI_COMM_WORLD);
			MPI_Send(proc_A->mat, (proc_A->ld_row)*(proc_A->ld_col), MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
		}
	}
}


void colToBuf(proc_mat *proc_A, double *col_block, int block, int k_proc)
{
	int block_length = proc_A->ld_row * block;
	int i;
	int j = 0;
	for(i = k_proc * block_length; i < (k_proc + 1)*block_length; ++i){
		col_block[j++] = proc_A->mat[i];
	}
}

void rowToBuf(proc_mat *proc_B, double *row_block, int block, int k_proc)
{
	int block_length = proc_B->ld_col*block;
	int i, j;
	i = k_proc * block;
	for(j = 0; j < block_length; ++j){
		row_block[j] = proc_B->mat[i];
		if((i+1) % block == 0)
			i = i + proc_B->ld_col-block + 1;
        else
			i++;
	}
}

void collect(int id, proc_mat *proc_C, proc_mat *old_C, int block)
{
	int global_i, global_j, i, j;
	global_i = (id/Q) * block;
	for(i = 0; i < proc_C->row; ++i){
		global_j = (id%Q) * block;
		for(j = 0; j < proc_C->col; ++j){
			old_C->mat[global_j * old_C->row+global_i] = proc_C->mat[j * proc_C->ld_row + i];
			if((j+1)%block == 0)
				global_j = global_j + block*(Q-1)+1;
			else
				global_j++;
		}
		if((i+1)%block == 0)
			global_i = global_i+block*(P-1)+1;
		else
			global_i++;
	}
}

void* mmm(void* param){
	thread_data* data = (thread_data*)param;
	data->thread_B = data->thread_B+(data->id)*(block*data->thread_col);
	data->thread_C = data->thread_C+(data->id)*(data->thread_row)*(data->thread_col);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, data->thread_row, data->thread_col, block, 1, data->thread_A, data->thread_row, data->thread_B, block, 1, data->thread_C, data->thread_row);
    return NULL;
}

void pdgemm(double *A, double *B, double *C, int n){

	proc_mat *proc_A, *proc_B, *proc_C, *old_A, *old_B, *old_C;
	MPI_Status status;

	int nb;
	nb = (n+block-1)/block;

	proc_A = (proc_mat*)malloc(sizeof(proc_mat));
	proc_A->row = 0;
	proc_A->col = 0;
	proc_A->b_row = (nb+P-1)/P;
	proc_A->b_col = (nb+Q-1)/Q;
	proc_A->ld_row = proc_A->b_row * block;
	proc_A->ld_col = proc_A->b_col * block;
	double *local_A = (double*)malloc(sizeof(double)*(proc_A->ld_row)*(proc_A->ld_col));
	proc_A->mat = local_A;

	proc_B = (proc_mat*)malloc(sizeof(proc_mat));
	proc_B->row = 0;
	proc_B->col = 0;
	proc_B->b_row = (nb+P-1)/P;
	proc_B->b_col = (nb+Q-1)/Q;
	proc_B->ld_row = proc_B->b_row * block;
	proc_B->ld_col = proc_B->b_col * block;
	double *local_B = (double*)malloc(sizeof(double)*(proc_B->ld_row)*(proc_B->ld_col));
	proc_B->mat = local_B;
	
	proc_C = (proc_mat*)malloc(sizeof(proc_mat));
	proc_C->row = 0;
	proc_C->col = 0;
	proc_C->b_row = proc_A->b_row;
	proc_C->b_col = proc_B->b_col;
	proc_C->ld_row = proc_A->ld_row;
	proc_C->ld_col = proc_B->ld_col;
	double *local_C = (double*)malloc(sizeof(double)*(proc_C->ld_row)*(proc_C->ld_col));
	proc_C->mat = local_C;

	MPI_Barrier(MPI_COMM_WORLD);		

	if(rank == 0){
		old_A = (proc_mat*)malloc(sizeof(proc_mat));
		old_A->row = n;
		old_A->col = n;
		old_A->b_row = nb;
		old_A->b_col = nb;
		old_A->ld_row = n;
		old_A->ld_col = n;
		old_A->mat = A;

		old_B = (proc_mat*)malloc(sizeof(proc_mat));
		old_B->row = n;
		old_B->col = n;
		old_B->b_row = nb;
		old_B->b_col = nb;
		old_B->ld_row = n;
		old_B->ld_col = n;
		old_B->mat = B;

		old_C = (proc_mat*)malloc(sizeof(proc_mat));
		old_C->row = n;
		old_C->col = n;
		old_C->b_row = nb;
		old_C->b_col = nb;
		old_C->ld_row = n;
		old_C->ld_col = n;
		old_C->mat = C;


		// matrix be distributed on differnet processes
		cyclicData(proc_A, old_A);
		cyclicData(proc_B, old_B);		
		cyclicData(proc_C, old_C);

	}
    else
    {
		// 1 to size-1 receive
		// tag = 0
		MPI_Recv(proc_A, sizeof(proc_mat), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(local_A, (proc_A->ld_row)*(proc_A->ld_col), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		proc_A->mat = local_A;
				
		MPI_Recv(proc_B, sizeof(proc_mat), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(local_B, (proc_B->ld_row)*(proc_B->ld_col), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		proc_B->mat = local_B;
		
		MPI_Recv(proc_C, sizeof(proc_mat), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(local_C, (proc_C->ld_row)*(proc_C->ld_col), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		proc_C->mat = local_C;
	}	

    MPI_Barrier(MPI_COMM_WORLD);			

	proc_C->row = proc_A->row;
	proc_C->col = proc_B->col;
	
	double* row_block=(double*)malloc((proc_B->ld_col)*block*sizeof(double));
	double* col_block=(double*)malloc((proc_A->ld_row)*block*sizeof(double));

	/* Split comm into row and column comms */ 
	MPI_Comm row_world, col_world;
	/* color by row, rank by column */
	MPI_Comm_split(MPI_COMM_WORLD, rank/Q, rank%Q, &row_world);
	/* color by column, rank by row */ 
	MPI_Comm_split(MPI_COMM_WORLD, rank%Q, rank/Q, &col_world);
	
	// SUMMA
	int id_row;				// rank in row comm
	int id_col;				// rank in column comm

	MPI_Comm_rank(row_world, &id_row);
	MPI_Comm_rank(col_world, &id_col);

	int k_index;			// block index of k on whole matrix
	int k_proc;				// block index of a local process: from 0 to nb-1
	int send_p;				// current block k on which process

	for(k_index = 0; k_index < nb; ++k_index){
		// broadcast A( , k), kth column to other processes
		send_p = k_index%Q;
		k_proc = k_index/Q + k_index%Q - id_row;
		if(id_row == send_p){
			colToBuf(proc_A, col_block, block, k_proc);	
		}
		MPI_Bcast(col_block, proc_A->ld_row*block, MPI_DOUBLE, send_p, row_world);
		
		// broadcast B(k,:), kth row to other processes
		send_p = k_index%P;
		k_proc = k_index/P + k_index%P - id_col;
		if(id_col == send_p){
			rowToBuf(proc_B, row_block, block, k_proc);	
		}
		MPI_Bcast(row_block, proc_B->ld_col*block, MPI_DOUBLE, send_p, col_world);
		
	    // do the dgemm using pthread
		pthread_t tid[NUM_THREADS];
		int j;
		for(j = 0; j < NUM_THREADS; ++j)
        {
			thread_data* data = (thread_data*)malloc(sizeof(thread_data));
			data->id = j;
			data->thread_row = proc_A->ld_row;
			data->thread_col = proc_B->ld_col/NUM_THREADS;
			data->thread_A = col_block;
			data->thread_B = row_block;
			data->thread_C = local_C;
			pthread_create(&(tid[j]), NULL, mmm, (void *)data);
		}
		for(j=0; j < NUM_THREADS; ++j){
			pthread_join(tid[j], NULL);
		}
	}

	// Collect result from other processes to process 0
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
    {
		proc_C->mat = local_C;
		collect(0, proc_C, old_C, block);
		int id;
		for(id = 1; id < size; ++id){
			// tag = 0
			MPI_Recv(proc_C, sizeof(proc_mat), MPI_BYTE, id, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(local_C, proc_C->ld_row*proc_C->ld_col, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, &status);
			proc_C->mat = local_C;
			collect(id, proc_C, old_C, block);
		}	
	}
    else
    {
		MPI_Send(proc_C, sizeof(proc_mat), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(local_C, proc_C->ld_row * proc_C->ld_col, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	free(proc_A);
	free(proc_B);
	free(proc_C);
}
