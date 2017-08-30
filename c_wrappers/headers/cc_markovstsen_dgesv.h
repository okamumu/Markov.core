#ifndef _CC_MARKOVSTSEN_DGESV_H
#define _CC_MARKOVSTSEN_DGESV_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_markovstsen_dgesv_dense(int n, double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* s, int incs, int* info);

void cc_markovstsen_dgesv_csr(int n,
  double* spQ, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info);

void cc_markovstsen_dgesv_csc(int n,
  double* spQ, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info);

void cc_markovstsen_dgesv_coo(int n,
  double* spQ, int* rowind0, int* colind0, int nnz0,
  double* dQ, int* rowind1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info);

#ifdef __cplusplus
}
#endif
#endif
