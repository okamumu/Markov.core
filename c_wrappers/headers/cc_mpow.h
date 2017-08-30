#ifndef _CC_MPOW_H
#define _CC_MPOW_H

#ifdef __cplusplus
extern "C" {
#endif

////

void cc_mpow_dense(char trans, int n, double* MA, int lda, double* ME, int lde, int m, int* info);

void cc_mpow_csr(char trans, int n, double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, int m, int* info);

void cc_mpow_csc(char trans, int n, double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, int m, int* info);

void cc_mpow_coo(char trans, int n, double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, int m, int* info);

#ifdef __cplusplus
}
#endif
#endif
