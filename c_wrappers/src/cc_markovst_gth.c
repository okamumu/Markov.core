///

void f90_markovst_gth_dense_(int* n, double* Q, int* ldq, double* x, int* incx);

void f90_markovst_gth_csr_(int* n, double* spQ, int* rowptr, int* colind, int* nnz,
  double* x, int* incx);

void f90_markovst_gth_csc_(int* n, double* spQ, int* colptr, int* rowind, int* nnz,
  double* x, int* incx);

void f90_markovst_gth_coo_(int* n, double* spQ, int* rowind, int* colind, int* nnz,
  double* x, int* incx);

///

void cc_markovst_gth_dense(int n, double* Q, int ldq, double* x, int incx) {

  f90_markovst_gth_dense_(&n, Q, &ldq, x, &incx);
}

void cc_markovst_gth_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_csr_(&n, spQ, rowptr, colind, &nnz, x, &incx);
}

void cc_markovst_gth_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_csc_(&n, spQ, colptr, rowind, &nnz, x, &incx);
}

void cc_markovst_gth_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_coo_(&n, spQ, rowind, colind, &nnz, x, &incx);
}
