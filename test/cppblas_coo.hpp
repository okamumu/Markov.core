/*
  blas.hpp
*/

namespace marlib {

  // basic functions: We expect VectorT = NumericVector and MatrixT = NumericMatrix in Rcpp

  // nnz
  template <>
  inline
  int dnnz(const coo_matrix& m) {
    return m.nnz();
  }

  // dcopy
  template <>
  inline
  coo_matrix& dcopy(const dense_matrix& x, coo_matrix& y) {
    dblas::dense_to_coo(nrow(x), ncol(x), &x[0], ld(x),
    &y[0], &y.rowind[0], &y.colind[0], y.nnz(), y.origin());
    return y;
  }

  template <>
  inline
  dense_matrix& dcopy(const coo_matrix& x, dense_matrix& y) {
    dblas::coo_to_dense(nrow(x), ncol(x), &x[0], &x.rowind[0], &x.colind[0], x.nnz(), x.origin(),
    &y[0], ld(y));
    return y;
  }

  // dgemv

  template <typename VectorT>
  inline
  VectorT& dgemvN(double alpha, const coo_matrix& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dcoomvN(nrow(A), ncol(A), alpha, &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  template <typename VectorT>
  inline
  VectorT& dgemvT(double alpha, const coo_matrix& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dcoomvT(nrow(A), ncol(A), alpha, &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  // BLAS level 3

  // dgemm

  template <typename MatrixT>
  inline
  MatrixT& dgemmNN(double alpha, const coo_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcoommNN(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmNT(double alpha, const coo_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcoommNT(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmTN(double alpha, const coo_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcoommTN(nrow(C), ncol(C), nrow(A), alpha,  &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmTT(double alpha, const coo_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcoommTT(nrow(C), ncol(C), nrow(A), alpha,  &A[0], &A.rowind[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

}
