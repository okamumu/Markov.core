#ifndef _CC_MEXPCONV_UNIF_H
#define _CC_MEXPCONV_UNIF_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_mexpconv_unif_dense_vec(char transQ, char transH,
  int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* H, int ldh, double atol);

void cc_mexpconv_unif_csr_vec(char transQ, char transH,
  int n, double* spP, int* rowptr, int* colind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol);

void cc_mexpconv_unif_csc_vec(char transQ, char transH,
  int n, double* spP, int* colptr, int* rowind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol);

void cc_mexpconv_unif_coo_vec(char transQ, char transH,
  int n, double* spP, int* rowind, int* colind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol);

#ifdef __cplusplus
}
#endif
#endif
