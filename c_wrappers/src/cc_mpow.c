
void f90_mpow_dense_(char* trans, int* n, double* MA, int* lda, double* ME, int* lde, int* m, int* info);

void f90_mpow_csr_(char* trans, int* n, double* spMA, int* rowptr, int* colind, int* nnz,
  double* ME, int* lde, int* m, int* info);

void f90_mpow_csc_(char* trans, int* n, double* spMA, int* colptr, int* rowind, int* nnz,
  double* ME, int* lde, int* m, int* info);

void f90_mpow_coo_(char* trans, int* n, double* spMA, int* rowind, int* colind, int* nnz,
  double* ME, int* lde, int* m, int* info);

//

void cc_mpow_dense(char trans, int n, double* MA, int lda, double* ME, int lde, int m, int* info) {
  f90_mpow_dense_(&trans, &n, MA, &lda, ME, &lde, &m, info);
}

void cc_mpow_csr(char trans, int n, double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_csr_(&trans, &n, spMA, rowptr, colind, &nnz, ME, &lde, &m, info);
}

void cc_mpow_csc(char trans, int n, double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_csc_(&trans, &n, spMA, colptr, rowind, &nnz, ME, &lde, &m, info);
}

void cc_mpow_coo(char trans, int n, double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_coo_(&trans, &n, spMA, rowind, colind, &nnz, ME, &lde, &m, info);
}
