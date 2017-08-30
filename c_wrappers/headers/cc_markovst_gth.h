#ifndef _CC_MARKOVST_GTH_H
#define _CC_MARKOVST_GTH_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_markovst_gth_dense(int n, double* Q, int ldq, double* x, int incx);

void cc_markovst_gth_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int incx);

void cc_markovst_gth_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int incx);

void cc_markovst_gth_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* x, int incx);

#ifdef __cplusplus
}
#endif
#endif
