#ifndef _CC_UNIF_MATRIX_H
#define _CC_UNIF_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_unif_dense(int n, double* Q, int ldq, double* P, int ldp, double* qv, double ufact);

void cc_unif_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* spP, double* qv, double ufact);

void cc_unif_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* spP, double* qv, double ufact);

void cc_unif_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* spP, double* qv, double ufact);

#ifdef __cplusplus
}
#endif
#endif
