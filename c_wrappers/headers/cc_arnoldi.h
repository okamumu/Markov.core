#ifndef _CC_ARNOLDI_H
#define _CC_ARNOLDI_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_arnoldi_dense(char trans, int n, double* A, int lda,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info);

void cc_arnoldi_csr(char trans, int n, double* A, int* rowptr, int* colind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info);

void cc_arnoldi_csc(char trans, int n, double* A, int* colptr, int* rowind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info);

void cc_arnoldi_coo(char trans, int n, double* A, int* rowind, int* colind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info);

#ifdef __cplusplus
}
#endif
#endif
