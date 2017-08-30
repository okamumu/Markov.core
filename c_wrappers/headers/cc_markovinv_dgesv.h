#ifndef _CC_MARKOVINV_DGESV_H
#define _CC_MARKOVINV_DGESV_H

#ifdef __cplusplus
extern "C" {
#endif

////

void cc_markovinv_dgesv_dense(char trans, int n, int nrhs,
  double* Q, int ldq,
  double* x, int ldx, double* y, int ldy, int* info);

void cc_markovinv_dgesv_csr(char trans, int n, int nrhs,
  double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info);

void cc_markovinv_dgesv_csc(char trans, int n, int nrhs,
  double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info);

void cc_markovinv_dgesv_coo(char trans, int n, int nrhs,
  double* spQ, int* rowind, int* colind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info);

#ifdef __cplusplus
}
#endif
#endif
