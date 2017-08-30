void f90_markovst_dgesv_dense_(int* n, double* Q, int* ldq,
  double* x, int* incx, int* info);

void f90_markovst_dgesv_csr_(int* n, double* spQ, int* rowptr, int* colind, int* nnz,
  double* x, int* incx, int* info);

void f90_markovst_dgesv_csc_(int* n, double* spQ, int* colptr, int* rowind, int* nnz,
  double* x, int* incx, int* info);

void f90_markovst_dgesv_coo_(int* n, double* spQ, int* rowind, int* colind, int* nnz,
  double* x, int* incx, int* info);

//

void cc_markovst_dgesv_dense(int n, double* Q, int ldq,
  double* x, int incx, int* info) {

  f90_markovst_dgesv_dense_(&n, Q, &ldq, x, &incx, info);
}

void cc_markovst_dgesv_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int incx, int* info) {

  f90_markovst_dgesv_csr_(&n, spQ, rowptr, colind, &nnz, x, &incx, info);
}

void cc_markovst_dgesv_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int incx, int* info) {

  f90_markovst_dgesv_csc_(&n, spQ, colptr, rowind, &nnz, x, &incx, info);
}

void cc_markovst_dgesv_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* x, int incx, int* info) {

  f90_markovst_dgesv_coo_(&n, spQ, rowind, colind, &nnz, x, &incx, info);
}
