#ifndef _CC_MARKOVST_SOR_H
#define _CC_MARKOVST_SOR_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_markovst_sor_dense(int n, double* Q, int ldq, double* xstart, int incxs,
  double* x, int incx, int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*));

void cc_markovst_sor_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*));

void cc_markovst_sor_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*));

#ifdef __cplusplus
}
#endif
#endif
