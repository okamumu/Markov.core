/*
  blas.hpp
*/

namespace marlib {

  // basic functions: We expect VectorT = NumericVector and MatrixT = NumericMatrix in Rcpp

  // nnz
  template <>
  inline
  int dnnz(const csr_matrix& m) {
    return m.nnz();
  }

  // dcopy
  template <>
  inline
  csr_matrix& dcopy(const dense_matrix& x, csr_matrix& y) {
    dblas::dense_to_csr(nrow(x), ncol(x), &x[0], ld(x),
    &y[0], &y.rowptr[0], &y.colind[0], y.nnz(), y.origin());
    return y;
  }

  template <>
  inline
  dense_matrix& dcopy(const csr_matrix& x, dense_matrix& y) {
    dblas::csr_to_dense(nrow(x), ncol(x), &x[0], &x.rowptr[0], &x.colind[0], x.nnz(), x.origin(),
    &y[0], ld(y));
    return y;
  }

  // dgemv

  template <typename VectorT>
  inline
  VectorT& dgemvN(double alpha, const csr_matrix& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dcsrmvN(nrow(A), ncol(A), alpha, &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  template <typename VectorT>
  inline
  VectorT& dgemvT(double alpha, const csr_matrix& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dcsrmvT(nrow(A), ncol(A), alpha, &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  // BLAS level 3

  // dgemm

  template <typename MatrixT>
  inline
  MatrixT& dgemmNN(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmNN(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmNT(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmNT(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmTN(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmTN(nrow(C), ncol(C), nrow(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT>
  inline
  MatrixT& dgemmTT(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmTT(nrow(C), ncol(C), nrow(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

}
