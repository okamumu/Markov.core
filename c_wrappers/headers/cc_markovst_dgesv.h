#ifndef _CC_MARKOVST_DGESV_H
#define _CC_MARKOVST_DGESV_H

#ifdef __cplusplus
extern "C" {
#endif


void cc_markovst_dgesv_dense(int n, double* Q, int ldq,
  double* x, int incx, int* info);

void cc_markovst_dgesv_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int incx, int* info);

void cc_markovst_dgesv_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int incx, int* info);

void cc_markovst_dgesv_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* x, int incx, int* info);

#ifdef __cplusplus
}
#endif
#endif
