
void f90_markovqst_gs_dense_(int* n, double* Q, int* ldq, double* xi, int* incxi,
  double* xstart, int* incxs, double* x, int* incx,
  double* gam, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovqst_gs_csr_(int* n, double* Q, int* rowptr, int* colind, int* nnz,
  double* xi, int* incxi, double* xstart, int* incxs, double* x, int* incx,
  double* gam, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovqst_gs_csc_(int* n, double* Q, int* colptr, int* rowind, int* nnz,
  double* xi, int* incxi, double* xstart, int* incxs, double* x, int* incx,
  double* gam, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

//

void cc_markovqst_gs_dense(int n, double* Q, int ldq, double* xi, int incxi,
  double* xstart, int incxs, double* x, int incx,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst_gs_dense_(&n, Q, &ldq, xi, &incxi,
    xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovqst_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst_gs_csr_(&n, Q, rowptr, colind, &nnz, xi, &incxi,
    xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovqst_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst_gs_csc_(&n, Q, colptr, rowind, &nnz, xi, &incxi,
    xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}
