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

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTN(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmTN(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTT(double alpha, const csr_matrix& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dcsrmmTT(nrow(C), ncol(C), ncol(A), alpha,  &A[0], &A.rowptr[0], &A.colind[0], A.nnz(), A.origin(),
    &B[0], ld(B), beta, &C[0], ld(C));
    return C;
  }

  ///

  // template <typename MatrixT>
  // inline
  // Rcpp::NumericVector& dgemmNN(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
  //   double beta, Rcpp::NumericVector& C) {
  //   return dgemvN(alpha, A, B, beta, C);
  // }
  //
  // template <typename MatrixT>
  // inline
  // Rcpp::NumericVector& dgemmNT(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
  //   double beta, Rcpp::NumericVector& C) {
  //   return dgemvN(alpha, A, B, beta, C);
  // }
  //
  // template <typename MatrixT>
  // inline
  // Rcpp::NumericVector& dgemmTN(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
  //   double beta, Rcpp::NumericVector& C) {
  //   return dgemvT(alpha, A, B, beta, C);
  // }
  //
  // template <typename MatrixT>
  // inline
  // Rcpp::NumericVector& dgemmTT(double alpha, const MatrixT& A, const Rcpp::NumericVector& B,
  //   double beta, Rcpp::NumericVector& C) {
  //   return dgemvT(alpha, A, B, beta, C);
  // }

//   // lapack
//
//   // dgesv: solve A X = B
//   // The solution is assgined to B
//
//   template <typename ValueT>
//   dense_matrix<ValueT>& dgesv(
//     const dense_matrix<ValueT>& A, dense_matrix<ValueT>& B);
//
//   template <typename ValueT>
//   vector<ValueT>& dgesv(
//     const dense_matrix<ValueT>& A, vector<ValueT>& B);
//
// #ifdef F77BLAS
//   template <>
//   inline
//   dense_matrix<double>& dgesv(
//     const dense_matrix<double>& A, dense_matrix<double>& B) {
//     assert(A.nrow() == A.ncol());
//     dense_matrix<double> MA = A.clone();
//     array<int> ipiv(A.nrow());
//     int info = dblas::dgesv(A.nrow(), B.ncol(), MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.ld());
//     assert(info == 0);
//     return B;
//   }
//
//   template <>
//   inline
//   vector<double>& dgesv(
//     const dense_matrix<double>& A, vector<double>& B) {
//     assert(A.nrow() == A.ncol());
//     dense_matrix<double> MA = A.clone();
//     array<int> ipiv(A.nrow());
//     int info = dblas::dgesv(A.nrow(), 1, MA.ptr(), MA.ld(), &ipiv[0], B.ptr(), B.size());
//     assert(info == 0);
//     return B;
//   }
// #endif

}
