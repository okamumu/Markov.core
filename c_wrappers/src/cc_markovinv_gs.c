//////

void f90_markovinv_gs_dense_(char* trans, int* n, int* nrhs,
  double* Q, int* ldq, double* x, int* ldx,
  double* ystart, int* ldys, double* y, int* ldy,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovinv_gs_csr_(char* trans, int* n, int* nrhs,
  double* Q, int* rowptr, int* colind, int* nnz, double* x, int* ldx,
  double* ystart, int* ldys, double* y, int* ldy,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovinv_gs_csc_(char* trans, int* n, int* nrhs,
  double* Q, int* colptr, int* rowind, int* nnz, double* x, int* ldx,
  double* ystart, int* ldys, double* y, int* ldy,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

//////

void cc_markovinv_gs_dense(char trans, int n, int nrhs,
  double* Q, int ldq, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovinv_gs_dense_(&trans, &n, &nrhs,
    Q, &ldq, x, &ldx, ystart, &ldys, y, &ldy,
    &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovinv_gs_csr(char trans, int n, int nrhs,
  double* Q, int* rowptr, int* colind, int nnz, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovinv_gs_csr_(&trans, &n, &nrhs,
    Q, rowptr, colind, &nnz, x, &ldx, ystart, &ldys, y, &ldy,
    &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovinv_gs_csc(char trans, int n, int nrhs,
  double* Q, int* colptr, int* rowind, int nnz, double* x, int ldx,
  double* ystart, int ldys, double* y, int ldy,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovinv_gs_csc_(&trans, &n, &nrhs,
    Q, colptr, rowind, &nnz, x, &ldx, ystart, &ldys, y, &ldy,
    &maxiter, &rtol, &steps, iter, rerror, info, callback);
}
