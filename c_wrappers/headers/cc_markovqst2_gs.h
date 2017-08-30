#ifndef _CC_MARKOVQST2_GS_H
#define _CC_MARKOVQST2_GS_H

#ifdef __cplusplus
extern "C" {
#endif


void cc_markovqst2_gs_dense(int n, double* Q, int ldq, double* xi, int incxi,
  double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovqst2_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovqst2_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

#ifdef __cplusplus
}
#endif
#endif
