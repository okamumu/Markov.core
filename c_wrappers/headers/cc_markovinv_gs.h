#ifndef _CC_MARKOVINV_GS_H
#define _CC_MARKOVINV_GS_H

#ifdef __cplusplus
extern "C" {
#endif

//////

void cc_markovinv_gs_dense(char trans, int n, int nrhs,
  double* Q, int ldq, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovinv_gs_csr(char trans, int n, int nrhs,
  double* Q, int* rowptr, int* colind, int nnz, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovinv_gs_csc(char trans, int n, int nrhs,
  double* Q, int* colptr, int* rowind, int nnz, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

#ifdef __cplusplus
}
#endif
#endif
