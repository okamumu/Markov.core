
// prototype

void f90_arnoldi_dense_(char* trans, int* n, double* A, int* lda,
  double* x, int* incx, int* m, double* H, int* ldh, double* V, int* ldv,
  double* beta, double* rnorm, double* tol, int* ite, int* info);

void f90_arnoldi_csr_(char* trans, int* n, double* A, int* rowptr, int* colind, int* nnz,
  double* x, int* incx, int* m, double* H, int* ldh, double* V, int* ldv,
  double* beta, double* rnorm, double* tol, int* ite, int* info);

void f90_arnoldi_csc_(char* trans, int* n, double* A, int* colptr, int* rowind, int* nnz,
  double* x, int* incx, int* m, double* H, int* ldh, double* V, int* ldv,
  double* beta, double* rnorm, double* tol, int* ite, int* info);

void f90_arnoldi_coo_(char* trans, int* n, double* A, int* rowind, int* colind, int* nnz,
  double* x, int* incx, int* m, double* H, int* ldh, double* V, int* ldv,
  double* beta, double* rnorm, double* tol, int* ite, int* info);

//

void cc_arnoldi_dense(char trans, int n, double* A, int lda,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info) {

  f90_arnoldi_dense_(&trans, &n, A, &lda, x, &incx, &m, H, &ldh, V, &ldv,
    beta, rnorm, &tol, &ite, info);
}

void cc_arnoldi_csr(char trans, int n, double* A, int* rowptr, int* colind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info) {

  f90_arnoldi_csr_(&trans, &n, A, rowptr, colind, &nnz,
    x, &incx, &m, H, &ldh, V, &ldv, beta, rnorm, &tol, &ite, info);
}

void cc_arnoldi_csc(char trans, int n, double* A, int* colptr, int* rowind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info) {

  f90_arnoldi_csc_(&trans, &n, A, colptr, rowind, &nnz,
    x, &incx, &m, H, &ldh, V, &ldv, beta, rnorm, &tol, &ite, info);
}

void cc_arnoldi_coo(char trans, int n, double* A, int* rowind, int* colind, int nnz,
  double* x, int incx, int m, double* H, int ldh, double* V, int ldv,
  double* beta, double* rnorm, double tol, int ite, int* info) {

  f90_arnoldi_coo_(&trans, &n, A, rowind, colind, &nnz,
    x, &incx, &m, H, &ldh, V, &ldv, beta, rnorm, &tol, &ite, info);
}
