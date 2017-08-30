#ifndef _CC_MARKOVST_GS_H
#define _CC_MARKOVST_GS_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_markovst_gs_dense(int n, double* Q, int ldq, double* xstart, int incxs,
  double* x, int incx, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovst_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovst_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

#ifdef __cplusplus
}
#endif
#endif
