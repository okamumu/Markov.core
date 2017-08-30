///

void f90_markovst_sor_dense_(int* n, double* Q, int* ldq, double* xstart, int* incxs,
  double* x, int* incx, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*));

void f90_markovst_sor_csr_(int* n, double* Q, int* rowptr, int* colind, int* nnz,
  double* xstart, int* incxs, double* x, int* incx,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*));

void f90_markovst_sor_csc_(int* n, double* Q, int* colptr, int* rowind, int* nnz,
  double* xstart, int* incxs, double* x, int* incx,
  int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*));

///

void cc_markovst_sor_dense(int n, double* Q, int ldq, double* xstart, int incxs,
  double* x, int incx, int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*)) {

  f90_markovst_sor_dense_(&n, Q, &ldq, xstart, &incxs, x, &incx,
    &maxiter, &rtol, &steps, iter, rerror, omega, info, callback, update_omega);
}

void cc_markovst_sor_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*)) {

  f90_markovst_sor_csr_(&n, Q, rowptr, colind, &nnz, xstart, &incxs,
    x, &incx, &maxiter, &rtol, &steps, iter, rerror, omega, info, callback, update_omega);
}

void cc_markovst_sor_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xstart, int incxs, double* x, int incx,
  int maxiter, double rtol, int steps,
  int* iter, double* rerror, double* omega, int* info,
  void (*callback)(int*,int*,double*,double*,int*,double*,double*,int*),
  void (*update_omega)(int*,int*,double*,double*,double*,double*,int*)) {

  f90_markovst_sor_csc_(&n, Q, colptr, rowind, &nnz, xstart, &incxs,
    x, &incx, &maxiter, &rtol, &steps, iter, rerror, omega, info, callback, update_omega);
}
