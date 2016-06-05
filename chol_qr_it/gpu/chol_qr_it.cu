#include <stdio.h>
#include <cutil.h>
#include "cublas.h"
#include "f2c.h"
#include "clapack.h"
#include "fcntl.h"
#include "cblas.h"
/*=============================================================================

extern "C" int sgemm_(char *, char *, int *, int *, int *, float *, float *, 
                      int *, float *, int *, float *, float *, int *);
extern "C" int sgesvd_(char *, char *, int *, int *, float *, int *, float *, 
                       float *, int *, float *, int *, float *, int *, int *);
extern "C" void sgeqrf_(int*, int*, float*, int*, float*, float*, int*, int*);
extern "C" int scopy_(int *, float*, int *, float*, int *);
extern "C" int strmm_(char*, char *, char*, char *, int *, int *, float *, 
                      float *, int *, float *, int *);
extern "C" int strsm_(char *, char *, char *, char *, int *, int *, 
                      float *, float *, int *, float *, int *);
extern "C" int ssyrk_(char *, char *, int *, int *, float *, float *, 
                      int *, float *, float *, int *);
//=============================================================================
*/
//void chol_qr_it(int m, int n, float *A, int lda, float *R){
void chol_qr_it(integer m, integer n, real *A, integer lda, real *R) {
	integer i = 0, k, j, info, lwork = n*n, n2 = n*n, one = 1;
    real *G, *U, *S, *VT, *vt, *tau, *work;  
    real cn = 200.f, alpha = 1.f, zero = 0.f, mins, maxs;
    
    G    = (float*)malloc(n * n * 4);
    VT   = (float*)malloc(n * n * 4);
    U   = (float*)malloc(n * n * 4);
    S    = (float*)malloc(    n * 4);
    work = (float*)malloc(lwork * 4);
    tau  = (float*)malloc(    n * 4);  

    do {
      i++;

			// G= A'*A
      sgemm_("t", "n", &n, &n, &m, &alpha, A, &m, A, &m, &zero, G, &n);
      //ssyrk_("l", "t", &n, &m, &alpha, A, &m, &zero, G, &n);
      //for(j=0; j<n; j++)
      //  for(k=0; k<j; k++)
      //     G[j*n+k] = G[k*n+j];

			// svd
      sgesvd_("n", "a", &n, &n, G, &n, S, U, &n, VT, &n, work, &lwork, &info);

      mins = 100.f, maxs = 0.f;
      for(k=0; k<n; k++){
        S[k] = sqrt(S[k]);

	if (S[k] < mins)  mins = S[k];
	if (S[k] > maxs)  maxs = S[k];
      }

      for(k=0; k<n;k++){
        vt = VT + k*n;
        for(j=0; j<n; j++)
          vt[j]*=S[j];
      } 
      sgeqrf_(&n, &n, VT, &n, tau, work, &lwork, &info);

      if (i==1)
        scopy_(&n2, VT, &one, R, &one);
      else
        strmm_("l", "u", "n", "n", &n, &n, &alpha, VT, &n, R, &n); // VT contains q and r, R is R

      strsm_("r", "u", "n", "n", &m, &n, &alpha, VT, &n, A, &m);    // A is Q

      if (mins > 0.00001f) 
        cn = maxs/mins;

      fprintf(stderr, "\nIteration %d, cond num = %f \n", i, cn);
    } while (cn > 100.f);

    free(G);
    free(VT);
    free(S);
    free(work);
    free(tau);
}

//=============================================================================

//void chol_qr_it_GPU(int m, int n, float *d_A, int lda, float *G, float *R, 
  //               float *work, int lwork){
void chol_qr_it_GPU(integer m, integer n, real *d_A, integer lda, real *G, real *R, 
                 real *work, integer lwork){
		integer i = 0, k, j, info, n2 = n*n, one = 1;
    real *U, *S, *VT, *d_VT, *vt, *tau, *h_G;
    real cn = 200.f, alpha = 1.f, zero = 0.f, mins, maxs;

    h_G  = (float*)malloc(n * n * 4);
    VT   = (float*)malloc(n * n * 4);
    U   = (float*)malloc(n * n * 4);
    cublasAlloc(n * n * 4, 1, (void **)&d_VT);
    S    = (float*)malloc(    n * 4);
    tau  = (float*)malloc(    n * 4);

    do {
      i++;
      cublasSgemm( 't', 'n', n, n, m, alpha, d_A, m, d_A, m, zero, G, n);
      cublasGetVector(n*n*4, 1, G, 1, h_G, 1);

      sgesvd_("n", "a", &n, &n, h_G, &n, S, U, &n, VT, &n, work, &lwork, &info);

      mins = 100.f, maxs = 0.f;
      for(k=0; k<n; k++){
        S[k] = sqrt(S[k]);

        if (S[k] < mins)  mins = S[k];
        if (S[k] > maxs)  maxs = S[k];
      }

      for(k=0; k<n;k++){
        vt = VT + k*n;
        for(j=0; j<n; j++)
          vt[j]*=S[j];
      }
      sgeqrf_(&n, &n, VT, &n, tau, work, &lwork, &info);

      if (i==1)
        scopy_(&n2, VT, &one, R, &one);
      else
        strmm_("l", "u", "n", "n", &n, &n, &alpha, VT, &n, R, &n);

      cublasSetVector(n*n*4, 1, VT, 1, d_VT, 1);
      cublasStrsm('r', 'u', 'n', 'n', m, n, alpha, d_VT, n, d_A, m);

      if (mins > 0.00001f)
        cn = maxs/mins;
      
      fprintf(stderr, "\nIteration %d, cond num = %f \n", i, cn);
    } while (cn > 100.f);

    free(VT);
    free(S);
    free(tau);
}

//=============================================================================
