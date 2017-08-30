///

void f90_markovinv_dgesv_dense_(char* trans, int* n, int* nrhs,
  double* Q, int* ldq,
  double* x, int* ldx, double* y, int* ldy, int* info);

void f90_markovinv_dgesv_csr_(char* trans, int* n, int* nrhs,
  double* spQ, int* rowptr, int* colind, int* nnz,
  double* x, int* ldx, double* y, int* ldy, int* info);

void f90_markovinv_dgesv_csc_(char* trans, int* n, int* nrhs,
  double* spQ, int* colptr, int* rowind, int* nnz,
  double* x, int* ldx, double* y, int* ldy, int* info);

void f90_markovinv_dgesv_coo_(char* trans, int* n, int* nrhs,
  double* spQ, int* rowind, int* colind, int* nnz,
  double* x, int* ldx, double* y, int* ldy, int* info);

//

void cc_markovinv_dgesv_dense(char trans, int n, int nrhs,
  double* Q, int ldq,
  double* x, int ldx, double* y, int ldy, int* info) {

  f90_markovinv_dgesv_dense_(&trans, &n, &nrhs,
    Q, &ldq, x, &ldx, y, &ldy, info);
}

void cc_markovinv_dgesv_csr(char trans, int n, int nrhs,
  double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info) {

  f90_markovinv_dgesv_csr_(&trans, &n, &nrhs,
    spQ, rowptr, colind, &nnz,
    x, &ldx, y, &ldy, info);
}

void cc_markovinv_dgesv_csc(char trans, int n, int nrhs,
  double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info) {

  f90_markovinv_dgesv_csc_(&trans, &n, &nrhs,
    spQ, colptr, rowind, &nnz,
    x, &ldx, y, &ldy, info);
}

void cc_markovinv_dgesv_coo(char trans, int n, int nrhs,
  double* spQ, int* rowind, int* colind, int nnz,
  double* x, int ldx, double* y, int ldy, int* info) {

  f90_markovinv_dgesv_coo_(&trans, &n, &nrhs,
    spQ, rowind, colind, &nnz,
    x, &ldx, y, &ldy, info);
}
