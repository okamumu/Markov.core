///

void f90_mexpint_unif_dense_vec_(char* trans, int* n, double* P, int* ldp, double* qv,
  int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* cy, int* inccy, double* atol);

void f90_mexpint_unif_csr_vec_(char* trans,
  int* n, double* spP, int* rowptr, int* colind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* cy, int* inccy, double* atol);

void f90_mexpint_unif_csc_vec_(char* trans,
  int* n, double* spP, int* colptr, int* rowind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* cy, int* inccy, double* atol);

void f90_mexpint_unif_coo_vec_(char* trans,
  int* n, double* spP, int* rowind, int* colind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* cy, int* inccy, double* atol);

void f90_mexpint_unif_dense_mat_(char* trans,
  int* n, double* P, int* ldp, double* qv,
  int* left, int* right, double* poi, double* weight,
  int* m, double* x, int* ldx, double* y, int* ldy, double* cy, int* ldcy, double* atol);

void f90_mexpint_unif_csr_mat_(char* trans,
  int* n, double* spP, int* rowptr, int* colind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  int* m, double* x, int* ldx, double* y, int* ldy, double* cy, int* ldcy, double* atol);

void f90_mexpint_unif_csc_mat_(char* trans,
  int* n, double* spP, int* colptr, int* rowind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  int* m, double* x, int* ldx, double* y, int* ldy, double* cy, int* ldcy, double* atol);

void f90_mexpint_unif_coo_mat_(char* trans,
  int* n, double* spP, int* rowind, int* colind, int* nnz, double* qv,
  int* left, int* right, double* poi, double* weight,
  int* m, double* x, int* ldx, double* y, int* ldy, double* cy, int* ldcy, double* atol);

///

void cc_mexpint_unif_dense_vec(char trans, int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* cy, int inccy, double atol) {

  f90_mexpint_unif_dense_vec_(&trans, &n, P, &ldp, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, cy, &inccy, &atol);
}

void cc_mexpint_unif_csr_vec(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* cy, int inccy, double atol) {

  f90_mexpint_unif_csr_vec_(&trans, &n, spP, rowptr, colind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, cy, &inccy, &atol);
}

void cc_mexpint_unif_csc_vec(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* cy, int inccy, double atol) {

  f90_mexpint_unif_csc_vec_(&trans, &n, spP, colptr, rowind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, cy, &inccy, &atol);
}

void cc_mexpint_unif_coo_vec(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* cy, int inccy, double atol) {

  f90_mexpint_unif_coo_vec_(&trans, &n, spP, rowind, colind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, cy, &inccy, &atol);
}

void cc_mexpint_unif_dense_mat(char trans,
  int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double* cy, int ldcy, double atol) {

  f90_mexpint_unif_dense_mat_(&trans, &n, P, &ldp, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, cy, &ldcy, &atol);
}

void cc_mexpint_unif_csr_mat(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double* cy, int ldcy, double atol) {

  f90_mexpint_unif_csr_mat_(&trans, &n, spP, rowptr, colind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, cy, &ldcy, &atol);
}

void cc_mexpint_unif_csc_mat(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double* cy, int ldcy, double atol) {

  f90_mexpint_unif_csc_mat_(&trans, &n, spP, colptr, rowind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, cy, &ldcy, &atol);
}

void cc_mexpint_unif_coo_mat(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double* cy, int ldcy, double atol) {

  f90_mexpint_unif_coo_mat_(&trans, &n, spP, rowind, colind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, cy, &ldcy, &atol);
}
