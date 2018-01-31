/*
  blas.hpp
*/

namespace marlib {

  // dfill

  template <typename VectorT>
  inline
  VectorT& dfill(VectorT& x, double alpha) {
    for (auto xptr = x.begin(); xptr != x.end(); xptr++) {
      *xptr = alpha;
    }
    return x;
  }

  // inc
  template <typename VectorT>
  inline
  int inc(const VectorT& x) {
    return 1;
  }

  // ld
  template <typename MatrixT>
  inline
  int ld(const Matrix& x) {
    return x.nrow();
  }

  // nnz
  template <typename MatrixT>
  inline
  int nnz(const MatrixT& m) {
    return dblas::nnz(m.nrow(), m.ncol(), &m[0], lda(m));
  }

  // dcopy
  template <typename VectorT, typename VectorT2>
  inline
  VectorT2& dcopy(const VectorT& x, VectorT2& y) {
    dblas::dcopy(x.size(), &x[0], inc(x), &y[0], inc(y));
    return y;
  }

  // BLAS level 1

  // daxpy

  template <typename VectorT>
  inline
  VectorT& daxpy(double alpha, const VectorT& x, VectorT& y) {
    dblas::daxpy(x.size(), alpha, &x[0], inc(x), &y[0], inc(y));
    return y;
  }

  // dscal

  template <typename VectorT>
  inline
  VectorT& dscal(double alpha, VectorT& x) {
    dblas::dscal(x.size(), alpha, &x[0], inc(x));
    return x;
  }

  // dasum

  template <typename VectorT>
  inline
  double dasum(const VectorT& x) {
    return dblas::dasum(x.size(), &x[0], inc(x));
  }

  // BLAS level 2

  // dgemv

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvN(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dgemv('N', A.nrow(), A.ncol(), alpha, &A[0], ld(A), &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvT(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dgemv('T', A.nrow(), A.ncol(), alpha, &A[0], ld(A), &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  // BLAS level 3

  // dgemm

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemm(char transA, char transB, double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    if (transA == 'N' || transA == 'n') {
      if (transB == 'N' || transB == 'n') {
        return dgemmNN(alpha, A, B, beta, C);
      } else {
        return dgemmNT(alpha, A, B, beta, C);
      }
    } else {
      if (transB == 'N' || transB == 'n') {
        return dgemmTN(alpha, A, B, beta, C);
      } else {
        return dgemmTT(alpha, A, B, beta, C);
      }
    }
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmNN(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('N', 'N', C.nrow(), C.ncol(), A.ncol(), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmNT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('N', 'T', C.nrow(), C.ncol(), A.ncol(), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTN(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'N', C.nrow(), C.ncol(), A.nrow(), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'T', C.nrow(), C.ncol(), A.nrow(), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
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
