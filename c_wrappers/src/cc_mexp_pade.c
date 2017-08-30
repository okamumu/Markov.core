///

void f90_mexp_pade_dense_(char* trans, int* n, double* alpha,
  double* MA, int* lda, double* ME, int* lde, double* eps);

void f90_mexp_pade_csr_(char* trans, int* n, double* alpha,
  double* spMA, int* rowptr, int* colind, int* nnz,
  double* ME, int* lde, double* eps);

void f90_mexp_pade_csc_(char* trans, int* n, double* alpha,
  double* spMA, int* colptr, int* rowind, int* nnz,
  double* ME, int* lde, double* eps);

void f90_mexp_pade_coo_(char* trans, int* n, double* alpha,
  double* spMA, int* rowind, int* colind, int* nnz,
  double* ME, int* lde, double* eps);

///

void cc_mexp_pade_dense(char trans, int n, double alpha,
  double* MA, int lda, double* ME, int lde, double eps) {

  f90_mexp_pade_dense_(&trans, &n, &alpha, MA, &lda, ME, &lde, &eps);
}

void cc_mexp_pade_csr(char trans, int n, double alpha,
  double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_csr_(&trans, &n, &alpha, spMA, rowptr, colind, &nnz,
    ME, &lde, &eps);
}

void cc_mexp_pade_csc(char trans, int n, double alpha,
  double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_csc_(&trans, &n, &alpha, spMA, colptr, rowind, &nnz,
    ME, &lde, &eps);
}

void cc_mexp_pade_coo(char trans, int n, double alpha,
  double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_coo_(&trans, &n, &alpha, spMA, rowind, colind, &nnz,
    ME, &lde, &eps);
}
