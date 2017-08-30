//

void f90_markovst_gs_dense_(int* n, double* Q, int* ldq, double* xstart, int* incxs,
  double* x, int* incx, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovst_gs_csr_(int* n, double* Q, int* rowptr, int* colind, int* nnz,
  double* xstart, int* incxs, double* x, int* incx,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovst_gs_csc_(int* n, double* Q, int* colptr, int* rowind, int* nnz,
  double* xstart, int* incxs, double* x, int* incx,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

///

void cc_markovst_gs_dense(int n, double* Q, int ldq, double* xstart, int incxs,
  double* x, int incx, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovst_gs_dense_(&n, Q, &ldq, xstart, &incxs, x, &incx,
    &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovst_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovst_gs_csr_(&n, Q, rowptr, colind, &nnz, xstart, &incxs,
    x, &incx, &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovst_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovst_gs_csc_(&n, Q, colptr, rowind, &nnz, xstart, &incxs,
    x, &incx, &maxiter, &rtol, &steps, iter, rerror, info, callback);
}
