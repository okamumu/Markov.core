#ifndef _CC_MEXP_UNIF_H
#define _CC_MEXP_UNIF_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_mexp_unif_dense_vec(char trans, int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol);

void cc_mexp_unif_csr_vec(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol);

void cc_mexp_unif_csc_vec(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol);

void cc_mexp_unif_coo_vec(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol);

void cc_mexp_unif_dense_mat(char trans, int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol);

void cc_mexp_unif_csr_mat(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol);

void cc_mexp_unif_csc_mat(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol);

void cc_mexp_unif_coo_mat(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol);

#ifdef __cplusplus
}
#endif
#endif
