/*
  blas.hpp
*/

namespace marlib {

  // basic functions: We expect VectorT = NumericVector and MatrixT = NumericMatrix in Rcpp

  // size
  template <typename VectorT>
  inline
  int size(const VectorT& x) {
    return x.size();
  }

  // inc
  template <typename VectorT>
  inline
  int inc(const VectorT& x) {
    return 1;
  }

  // nrow
  template <typename MatrixT>
  inline
  int nrow(const MatrixT& x) {
    return x.nrow();
  }

  // ncol
  template <typename MatrixT>
  inline
  int ncol(const MatrixT& x) {
    return x.ncol();
  }

  // ld
  template <typename MatrixT>
  inline
  int ld(const MatrixT& x) {
    return x.nrow();
  }

  // dfill

  template <typename VectorT>
  inline
  VectorT& dfill(VectorT& x, double alpha) {
    dblas::dfill(size(x), &x[0], 1, alpha);
    return x;
  }

  // nnz
  template <typename MatrixT>
  inline
  int dnnz(const MatrixT& m) {
    return dblas::dnnz(nrow(m), ncol(m), &m[0], ld(m));
  }

  // dcopy
  template <typename VectorT, typename VectorT2>
  inline
  VectorT2& dcopy(const VectorT& x, VectorT2& y) {
    dblas::dcopy(size(x), &x[0], inc(x), &y[0], inc(y));
    return y;
  }

  // BLAS level 1

  // daxpy

  template <typename VectorT>
  inline
  VectorT& daxpy(double alpha, const VectorT& x, VectorT& y) {
    dblas::daxpy(size(x), alpha, &x[0], inc(x), &y[0], inc(y));
    return y;
  }

  // dscal

  template <typename VectorT>
  inline
  VectorT& dscal(double alpha, VectorT& x) {
    dblas::dscal(size(x), alpha, &x[0], inc(x));
    return x;
  }

  // ddot

  template <typename VectorT, typename VectorT2>
  inline
  double dcopy(const VectorT& x, const VectorT2& y) {
    return dblas::ddot(size(x), &x[0], inc(x), &y[0], inc(y));
  }

  // dasum

  template <typename VectorT>
  inline
  double dasum(const VectorT& x) {
    return dblas::dasum(size(x), &x[0], inc(x));
  }

  // dnrm2

  template <typename VectorT>
  inline
  double dnrm2(const VectorT& x) {
    return dblas::dnrm2(size(x), &x[0], inc(x));
  }

  // damax

  template <typename VectorT>
  inline
  double damax(const VectorT& x) {
    return dblas::damax(size(x), &x[0], inc(x));
  }

  // BLAS level 2

  // dgemv

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvN(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dgemv('N', nrow(A), ncol(A), alpha, &A[0], ld(A), &x[0], inc(x), beta, &y[0], inc(y));
    return y;
  }

  template <typename VectorT, typename MatrixT>
  inline
  VectorT& dgemvT(double alpha, const MatrixT& A, const VectorT& x, double beta, VectorT& y) {
    dblas::dgemv('T', nrow(A), ncol(A), alpha, &A[0], ld(A), &x[0], inc(x), beta, &y[0], inc(y));
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
    dblas::dgemm('N', 'N', nrow(C), ncol(C), ncol(A), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmNT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('N', 'T', nrow(C), ncol(C), ncol(A), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTN(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'N', nrow(C), ncol(C), nrow(A), alpha, &A[0], ld(A), &B[0], ld(B),
      beta, &C[0], ld(C));
    return C;
  }

  template <typename MatrixT, typename MatrixT2>
  inline
  MatrixT& dgemmTT(double alpha, const MatrixT2& A, const MatrixT& B, double beta, MatrixT& C) {
    dblas::dgemm('T', 'T', nrow(C), ncol(C), nrow(A), alpha, &A[0], ld(A), &B[0], ld(B),
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
