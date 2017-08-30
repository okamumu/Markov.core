#ifndef _CC_BLAS_H
#define _CC_BLAS_H

#ifdef __cplusplus
extern "C" {
#endif


double f77blas_ddot(int n, double* x, int incx, double* y, int incy);

void f77blas_dcopy(int n, double* x, int incx, double* y, int incy);

void f77blas_daxpy(int n, double a, double* x, int incx, double* y, int incy);

#ifdef __cplusplus
}
#endif
#endif
