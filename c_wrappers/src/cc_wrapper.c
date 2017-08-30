
#include "f90_wrapper.h"

double f77blas_ddot(int n, double* x, int incx, double* y, int incy) {
  return ddot_(&n, x, &incx, y, &incy);
}

void f77blas_dcopy(int n, double* x, int incx, double* y, int incy) {
  dcopy_(&n, x, &incx, y, &incy);
}

void f77blas_daxpy(int n, double a, double* x, int incx, double* y, int incy) {
  daxpy_(&n, &a, x, &incx, y, &incy);
}

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

//////

void cc_markovqst2_gs_dense(int n, double* Q, int ldq, double* xi, int incxi,
  double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst2_gs_dense_(&n, Q, &ldq, xi, &incxi,
    xstart, &incxs, x, &incx, ystart, &incys, y, &incy,
    gam, &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovqst2_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst2_gs_csr_(&n, Q, rowptr, colind, &nnz, xi, &incxi,
    xstart, &incxs, x, &incx, ystart, &incys, y, &incy,
    gam, &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

void cc_markovqst2_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
  double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
  double* ystart, int incys, double* y, int incy,
  double* gam, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovqst2_gs_csc_(&n, Q, colptr, rowind, &nnz, xi, &incxi,
    xstart, &incxs, x, &incx, ystart, &incys, y, &incy,
    gam, &maxiter, &rtol, &steps, iter, rerror, info, callback);
}

///

// void cc_markovqst_gs_dense(int n, double* Q, int ldq, double* xi, int incxi,
//   double* xstart, int incxs, double* x, int incx,
//   double* gam, int maxiter, double rtol, int steps,
//   int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

//   f90_markovqst_gs_dense_(&n, Q, &ldq, xi, &incxi,
//     xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
//     iter, rerror, info, callback);
// }

// void cc_markovqst_gs_csr(int n, double* Q, int* rowptr, int* colind, int nnz,
//   double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
//   double* gam, int maxiter, double rtol, int steps,
//   int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

//   f90_markovqst_gs_csr_(&n, Q, rowptr, colind, &nnz, xi, &incxi,
//     xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
//     iter, rerror, info, callback);
// }

// void cc_markovqst_gs_csc(int n, double* Q, int* colptr, int* rowind, int nnz,
//   double* xi, int incxi, double* xstart, int incxs, double* x, int incx,
//   double* gam, int maxiter, double rtol, int steps,
//   int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

//   f90_markovqst_gs_csc_(&n, Q, colptr, rowind, &nnz, xi, &incxi,
//     xstart, &incxs, x, &incx, gam, &maxiter, &rtol, &steps,
//     iter, rerror, info, callback);
// }

///

// void cc_markovst_dgesv_dense(int n, double* Q, int ldq,
//   double* x, int incx, int* info) {

//   f90_markovst_dgesv_dense_(&n, Q, &ldq, x, &incx, info);
// }

// void cc_markovst_dgesv_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
//   double* x, int incx, int* info) {

//   f90_markovst_dgesv_csr_(&n, spQ, rowptr, colind, &nnz, x, &incx, info);
// }

// void cc_markovst_dgesv_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
//   double* x, int incx, int* info) {

//   f90_markovst_dgesv_csc_(&n, spQ, colptr, rowind, &nnz, x, &incx, info);
// }

// void cc_markovst_dgesv_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
//   double* x, int incx, int* info) {

//   f90_markovst_dgesv_coo_(&n, spQ, rowind, colind, &nnz, x, &incx, info);
// }

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

///

void cc_markovst_gth_dense(int n, double* Q, int ldq, double* x, int incx) {

  f90_markovst_gth_dense_(&n, Q, &ldq, x, &incx);
}

void cc_markovst_gth_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_csr_(&n, spQ, rowptr, colind, &nnz, x, &incx);
}

void cc_markovst_gth_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_csc_(&n, spQ, colptr, rowind, &nnz, x, &incx);
}

void cc_markovst_gth_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* x, int incx) {

  f90_markovst_gth_coo_(&n, spQ, rowind, colind, &nnz, x, &incx);
}

///

void cc_markovstsen_dgesv_dense(int n, double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_dense_(&n, Q, &ldq, dQ, &lddq, pis, &incp,
    s, &incs, info);
}

void cc_markovstsen_dgesv_csr(int n,
  double* spQ, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_csr_(&n, spQ, rowptr0, colind0, &nnz0,
    dQ, rowptr1, colind1, &nnz1, pis, &incp, s, &incs, info);
}

void cc_markovstsen_dgesv_csc(int n,
  double* spQ, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_csc_(&n, spQ, colptr0, rowind0, &nnz0,
    dQ, colptr1, rowind1, &nnz1, pis, &incp, s, &incs, info);
}

void cc_markovstsen_dgesv_coo(int n,
  double* spQ, int* rowind0, int* colind0, int nnz0,
  double* dQ, int* rowind1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_coo_(&n, spQ, rowind0, colind0, &nnz0,
    dQ, rowind1, colind1, &nnz1, pis, &incp, s, &incs, info);
}

///

void cc_markovstsen_gs_dense(int n,
  double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_dense_(&n, Q, &ldq, dQ, &lddq, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovstsen_gs_csr(int n,
  double* Q, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_csr_(&n, Q, rowptr0, colind0, &nnz0,
    dQ, rowptr1, colind1, &nnz1, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovstsen_gs_csc(int n,
  double* Q, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_csc_(&n, Q, colptr0, rowind0, &nnz0,
    dQ, colptr1, rowind1, &nnz1, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

///

void cc_mexp_pade_dense(char trans, int n, double alpha,
  double* MA, int lda, double* ME, int lde, double eps) {

  f90_mexp_pade_dense_(&trans, &n, &alpha, MA, &lda, ME, &lde, &eps);
}

void cc_mexp_pade_csr(char trans, int n, double alpha,
  double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_csr_(&trans, &n, &alpha, spMA, rowptr, colind, &nnz,
    ME, &lde, &eps);
}

void cc_mexp_pade_csc(char trans, int n, double alpha,
  double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_csc_(&trans, &n, &alpha, spMA, colptr, rowind, &nnz,
    ME, &lde, &eps);
}

void cc_mexp_pade_coo(char trans, int n, double alpha,
  double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, double eps) {

  f90_mexp_pade_coo_(&trans, &n, &alpha, spMA, rowind, colind, &nnz,
    ME, &lde, &eps);
}

///

void cc_mexp_unif_dense_vec(char trans, int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol) {

  f90_mexp_unif_dense_vec_(&trans, &n, P, &ldp, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, &atol);
}

void cc_mexp_unif_csr_vec(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol) {

  f90_mexp_unif_csr_vec_(&trans, &n, spP, rowptr, colind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, &atol);
}

void cc_mexp_unif_csc_vec(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol) {

  f90_mexp_unif_csc_vec_(&trans, &n, spP, colptr, rowind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, &atol);
}

void cc_mexp_unif_coo_vec(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double atol) {

  f90_mexp_unif_coo_vec_(&trans, &n, spP, rowind, colind, &nnz, &qv,
    &left, &right, poi, &weight, x, &incx, y, &incy, &atol);
}

void cc_mexp_unif_dense_mat(char trans, int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol) {

  f90_mexp_unif_dense_mat_(&trans, &n, P, &ldp, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, &atol);
}

void cc_mexp_unif_csr_mat(char trans,
  int n, double* spP, int* rowptr, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol) {

  f90_mexp_unif_csr_mat_(&trans, &n, spP, rowptr, colind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, &atol);
}

void cc_mexp_unif_csc_mat(char trans,
  int n, double* spP, int* colptr, int* rowind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol) {

  f90_mexp_unif_csc_mat_(&trans, &n, spP, colptr, rowind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, &atol);
}

void cc_mexp_unif_coo_mat(char trans,
  int n, double* spP, int* rowind, int* colind, int nnz, double qv,
  int left, int right, double* poi, double weight,
  int m, double* x, int ldx, double* y, int ldy, double atol) {

  f90_mexp_unif_coo_mat_(&trans, &n, spP, rowind, colind, &nnz, &qv,
    &left, &right, poi, &weight, &m, x, &ldx, y, &ldy, &atol);
}

///

void cc_mexpconv_unif_dense_vec(char transQ, char transH,
  int n, double* P, int ldp, double qv,
  int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* H, int ldh, double atol) {

  f90_mexpconv_unif_dense_vec_(&transQ, &transH, &n, P, &ldp, &qv,
    &left, &right, poi, &weight,
    x, &incx, y, &incy, z, &incz, H, &ldh, &atol);
}

void cc_mexpconv_unif_csr_vec(char transQ, char transH,
  int n, double* spP, int* rowptr, int* colind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol) {

  f90_mexpconv_unif_csr_vec_(&transQ, &transH, &n, spP, rowptr, colind, &nnz,
    &qv, &left, &right, poi, &weight, x, &incx, y, &incy, z, &incz, spH, &atol);
}

void cc_mexpconv_unif_csc_vec(char transQ, char transH,
  int n, double* spP, int* colptr, int* rowind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol) {

  f90_mexpconv_unif_csc_vec_(&transQ, &transH, &n, spP, colptr, rowind, &nnz,
    &qv, &left, &right, poi, &weight, x, &incx, y, &incy, z, &incz, spH, &atol);
}

void cc_mexpconv_unif_coo_vec(char transQ, char transH,
  int n, double* spP, int* rowind, int* colind, int nnz,
  double qv, int left, int right, double* poi, double weight,
  double* x, int incx, double* y, int incy, double* z, int incz,
  double* spH, int atol) {

  f90_mexpconv_unif_coo_vec_(&transQ, &transH, &n, spP, rowind, colind, &nnz,
    &qv, &left, &right, poi, &weight, x, &incx, y, &incy, z, &incz, spH, &atol);
}

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

///

int cc_poisson_rightbound(double lambda, double eps) {
  int right;
  f90_poisson_rightbound_(&lambda, &eps, &right);
  return right;
}

double cc_poisson_prob(double lambda, int left, int right, double* prob) {
  double weight;
  f90_poisson_prob_(&lambda, &left, &right, prob, &weight);
  return weight;
}

///

void cc_unif_dense(int n, double* Q, int ldq, double* P, int ldp, double* qv, double ufact) {
  f90_unif_dense_(&n, Q, &ldq, P, &ldp, qv, &ufact);
}

void cc_unif_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_csr_(&n, spQ, rowptr, colind, &nnz, spP, qv, &ufact);
}

void cc_unif_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_csc_(&n, spQ, colptr, rowind, &nnz, spP, qv, &ufact);
}

void cc_unif_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_coo_(&n, spQ, rowind, colind, &nnz, spP, qv, &ufact);
}

////

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

///

void cc_mpow_dense(char trans, int n, double* MA, int lda, double* ME, int lde, int m, int* info) {
  f90_mpow_dense_(&trans, &n, MA, &lda, ME, &lde, &m, info);
}

void cc_mpow_csr(char trans, int n, double* spMA, int* rowptr, int* colind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_csr_(&trans, &n, spMA, rowptr, colind, &nnz, ME, &lde, &m, info);
}

void cc_mpow_csc(char trans, int n, double* spMA, int* colptr, int* rowind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_csc_(&trans, &n, spMA, colptr, rowind, &nnz, ME, &lde, &m, info);
}

void cc_mpow_coo(char trans, int n, double* spMA, int* rowind, int* colind, int nnz,
  double* ME, int lde, int m, int* info) {

  f90_mpow_coo_(&trans, &n, spMA, rowind, colind, &nnz, ME, &lde, &m, info);
}
