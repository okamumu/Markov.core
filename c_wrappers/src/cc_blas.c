/// F77 blas

double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
void daxpy_(int* n, double* a, double* x, int* incx, double* y, int* incy);

double f77blas_ddot(int n, double* x, int incx, double* y, int incy) {
  return ddot_(&n, x, &incx, y, &incy);
}

void f77blas_dcopy(int n, double* x, int incx, double* y, int incy) {
  dcopy_(&n, x, &incx, y, &incy);
}

void f77blas_daxpy(int n, double a, double* x, int incx, double* y, int incy) {
  daxpy_(&n, &a, x, &incx, y, &incy);
}
