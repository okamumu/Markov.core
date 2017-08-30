#ifndef _CC_MARKOVSTSEN_GS_H
#define _CC_MARKOVSTSEN_GS_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_markovstsen_gs_dense(int n,
  double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovstsen_gs_csr(int n,
  double* Q, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void cc_markovstsen_gs_csc(int n,
  double* Q, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

#ifdef __cplusplus
}
#endif
#endif
