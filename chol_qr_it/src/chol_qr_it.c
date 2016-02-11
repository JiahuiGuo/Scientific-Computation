/**********************************************
 * Author: Jiahui Guo
 * Email: guo.jiahui07@gmail.com
 * Department: EECS
 * CS594 Homework 3
 * Function: C version of Matlab code chol_qr_it
***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/f2c.h"
#include "../include/cblas.h"
#include "../include/clapack.h"

void printMatrix(double* matrix, int row, int col)
{
    printf("The matrix is:\n");
    size_t k;
    for (k = 0; k < row*col; ++k){
        printf("%f\t", matrix[k]);
    }
    printf("\n");
    size_t i, j;
    for (i = 0; i < row; ++i)
    {
        for (j = 0; j < col; ++j)
            printf("%f\t", matrix[i*row + j]);
        printf("\n");
    }
}


void printA(doublereal *A, int x, int y)
{
	int i, j;
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	for(i=0; i<x; i++)
	{
		for(j=0; j<y; j++)
		{
			printf("%f ", A[i*y+j]);
		}
		printf("\n");
	}
}

void majorConvert(doublereal *a,int n)
{
	doublereal temp;
	int i, j;
	for (i=0; i<n; i++)
		for(j=i; j<n; j++)
			{
				temp = a[i*n+j];
				a[i*n+j] = a[j*n+i];
				a[j*n+i] = temp;
			}
}

int main(int argc, char** argv)
{
	integer n, m;
	doublereal *x, *orgX;
	doublereal *u, *s, *v, *vt, *sigma;
	doublereal *q, *r;
	doublereal *G, *work, *work2;
	doublereal *RR, *RRtemp;
	doublereal cn =200;
	int i,j;
	doublereal workopt;
	integer lwork=-1, lwork2=-1;
	integer info;
	int index_max;

	doublereal *tempqr;
	doublereal *tau;
	n = 32;	//Rows
	m = 32; //Columns
	if (argc != 3)
	{
		printf("Please enter size of matrix like cholqr 6 4\n");
		exit(0);
	}
	else
	{
		m = atoi(argv[1]);
		n = atoi(argv[2]);
	}
	printf("m=%d n=%d\n",m,n);
    char all = 'A';
	char uplo = 'U';
	char diag = 'N';
	x = (doublereal *)malloc(m*n*sizeof(doublereal));
	orgX = (doublereal *)malloc(m*n*sizeof(doublereal));
	u = (doublereal *)malloc(n*n*sizeof(doublereal));
	s = (doublereal *)malloc(n*n*sizeof(doublereal));
	sigma = (doublereal *)malloc(n*sizeof(doublereal)); // smaller one between m and n
	v = (doublereal *)malloc(n*n*sizeof(doublereal));
	vt = (doublereal*)malloc(n*n*sizeof(doublereal));
	q = (doublereal *)malloc(m*n*sizeof(doublereal));
	r = (doublereal *)malloc(n*n*sizeof(doublereal));
	G = (doublereal *)malloc(n*n*sizeof(doublereal));
	tempqr = (doublereal *)malloc(n*n*sizeof(doublereal));
	tau = (doublereal *)malloc(n*sizeof(doublereal)); // smaller one between m and n
	RR = (doublereal *)malloc(n*n*sizeof(doublereal));
	RRtemp = (doublereal *)malloc(n*n*sizeof(doublereal));
	
	doublereal min, condition;
	for(i=0;i<m;i++)
		for(j=0; j<n;j++)
		{
			x[i*n+j] = i*n+j+1;
			orgX[i*n+j] = x[i*n+j];
		}
	// eye(RR)
	for (i=0; i<n; i++)
		for(j=0; j<n;j++)
		{
			if(i==j)
				RR[i*n+j] = 1;
			else 
				RR[i*n+j] = 0;
		}
	// 101 for CblasRowMajor, 112 for Trans X; alpha = 1; beta = 0
	cblas_dsyrk(101, 121, 112,n,m,1.0,x,n,0.0,G,n);

    while(cn > 100)
    {
	for (i=0; i<n; i++)
		for(j=i+1; j<n; j++)
			G[j*n+i] = G[i*n+j];
	lwork=-1;
	dgesvd_(&all, &all, &n, &n, G, &n, sigma, u, &n, vt, &n, &workopt, &lwork, &info);
	lwork=(integer)workopt;
	work = (doublereal *)malloc(lwork*sizeof(doublereal));
	dgesvd_(&all, &all, &n, &n, G, &n, sigma, u, &n, vt, &n, work, &lwork, &info);
	if(info >0 )
		{
			printf("SVD failed to converge!!\n");
			exit(1);
		}
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			if(i==j)
				s[i*n+j] = sqrt(sigma[i]);
			else
				s[i*n+j] = 0.0;
		}
	// [q,r]=qr(sqrt(s)*v')
	// RowMajor, sqrt(s) NoTrans, vt NoTrans
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1, s, n, vt, n, 0, tempqr, n);
	// [q,r]=qr(tempqr);
	work2 =(doublereal *)malloc(n*sizeof(doublereal));
	doublereal temp;
	// attention to the col/row major issues!!
	//printA(tempqr, m,n);
	for(i=0; i<n; i++)
		for(j=i; j<n;j++)
		{
			temp = tempqr[i*n+j];
			tempqr[i*n+j] = tempqr[j*n+i];
			tempqr[j*n+i] = temp;
		}
	dgeqrf_(&n, &n, tempqr, &n, tau, work2, &n, &info);
	for (i=0; i<n; i++)
		for(j=0; j<n; j++)
		{
			if(i<=j)
				r[i*n+j] = tempqr[j*n+i];
			else
				r[i*n+j] = 0;
		}
	// +/- of r is different from Matlab result. 
	// cn = sqrt(cond(s))
	// cond(s)
	index_max = cblas_idamax(n,sigma,1);
	min = sigma[0];
	for (i=1;i<n; i++)
	{
	//	printf("sigma[%d] is %f\n", i, sigma[i]);
		if(sigma[i]<min)
			min = sigma[i];
	}
	condition = sqrt(sigma[index_max]/min);
	//printf("Condition = %f, max is %f min is %f\n", condition, sigma[index_max], min);
	cn = condition;
	// r = inv (r)
	// I have manually changed the +/- of r[1][*]; 
	// Check this later
	for (i=0;i<n;i++)
		r[n+i]=-r[n+i];
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n,n,n,1,r,n,RR,n,0,RRtemp,n);
	for (i=0; i<n; i++)
		for (j=0; j<n;j++)
			RR[i*n+j] = RRtemp[i*n+j];
	// rowmajor -> colmajor
	majorConvert(r, n);
	// inv(r)
	dtrtri_(&uplo, &diag, &n, r, &n, &info);
	// ColMajor to RowMajor
	majorConvert(r, n);
	// Q = Q*inv(r);-> x=x*r
	doublereal *Xtemp;
	Xtemp = (doublereal *)malloc(m*n*sizeof(doublereal));
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, n, 1, x, n, r, n, 0, Xtemp, n);

	for (i=0; i<m; i++)
		for(j=0; j<n; j++)
			x[i*n+j] = Xtemp[i*n+j];
	free(Xtemp);

	if (cn > 100)
		cblas_dsyrk(101, 121, 112,n,m,1.0,x,n,0.0,G,n);
		// Get x'*x
		// 101 for CblasRowMajor, 112 for Trans X; alpha = 1; beta = 0

	printf("***************************\n");
	printf("cn is %f\n", cn);
    }//end of while(cn>200)
	printf("***************************\n");
	printf("cn is %f\n", cn);
	printf("X is \n");
	printA(orgX, m,n);
	printf("Q is \n");
	printA(x, m, n);
	printf("R is \n");
	printA(RR, n, n);

	// I-Q'*Q
	doublereal *test;
	test = (doublereal *)malloc(n*n*sizeof(doublereal));
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			{
				if(i==j)
					test[i*n+j] = 1.0;
				else
					test[i*n+j] = 0.0;
			}
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n,n,m,1,x,n,x,n,-1,test,n);
	//printA(test,n,n);
	doublereal normQinvQ; 
	doublereal norm2;
	all = '1';
	normQinvQ = dlange_(&all, &n, &n, test, &n, work);
	doublereal	temp =0;
	for(i=0; i<n; i++)
		for(j=0;j<n; j++)
			temp = temp+test[i*n+j]	;
	printf("\nnorm Q'*Q (%d by %d) from dlange_ is %f\n", n,n, normQinvQ);
	printf("I-q'*q is %f\n", temp);
	free(test);
	// X-q*r
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m,n,n,1,x,n,RR,n,-1,orgX,n);
	norm2 = dlange_(&all, &m, &n, orgX, &m, work);
	temp=0;
	for(i=0;i<m;i++)
		for(j=0; j<n;j++)
			temp = temp + orgX[i*n+j];
	//printA(orgX, m,n);
	printf("Norm x-q*r (%d by %d) from dlange_ is %f\n", m,n,norm2);
	printf("x-q*r is %f\n", temp);

	free(work2);
	free(orgX);
	free(x);
	free(u);
	free(s);
	free(v);
	free(vt);
	free(q);
	free(r);
	free(G);
	free(RR);
	free(RRtemp);
    
    return EXIT_SUCCESS;
}
