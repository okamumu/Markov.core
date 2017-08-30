#ifndef _CC_MEXP_PADE_H
#define _CC_MEXP_PADE_H

#ifdef __cplusplus
extern "C" {
#endif

///

void cc_mexp_pade_dense(char trans, int n, double alpha,
  double* MA, int lda, double* ME, int lde, double eps);

void cc_mexp_pade_csr(char trans, int n, double alpha,
  double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, double eps);

void cc_mexp_pade_csc(char trans, int n, double alpha,
  double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, double eps);

void cc_mexp_pade_coo(char trans, int n, double alpha,
  double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, double eps);

#ifdef __cplusplus
}
#endif
#endif
