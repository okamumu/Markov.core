/*
  blas.h
  wrapper for blas
 */

namespace dblas {

  // level 2

  inline void dcsrmvN(const int m, const int n, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *x, const int incx, const double beta, double *y, const int incy) {
    dscal(m, beta, y, incy);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        y[i*incy] += alpha * A[z] * x[j*incx];
      }
    }
  }

  inline void dcsrmvT(const int m, const int n, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *x, const int incx, const double beta, double *y, const int incy) {
    dscal(n, beta, y, incy);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        y[j*incy] += alpha * A[z] * x[i*incx];
      }
    }
  }

  inline void dcsrr(const int m, const int n, const double alpha,
    const double *x, const int incx, const double *y, const int incy,
    double *A, const int *rowptr, const int *colind, const int nnz, const int origin) {
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        A[z] += alpha * x[i*incx] * y[j*incy];
      }
    }
  }

  inline void dcsrmmNN(const int m, const int n, const int k, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *B, const int ldb, const double beta, double *C, const int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = &B[0];
        double* Cptr = &C[0];
        for (int v=0; v<k; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j];
        }
      }
    }
  }

  inline void dcsrmmTN(const int m, const int n, const int k, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *B, const int ldb, const double beta, double *C, const int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<k; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i];
        }
      }
    }
  }

  inline void dcsrmmNT(const int m, const int n, const int k, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *B, const int ldb, const double beta, double *C, const int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = &B[0];
        double* Cptr = &C[0];
        for (int v=0; v<k; v++, Bptr+=1, Cptr+=1) {
          Cptr[i*ldc] += alpha * A[z] * Bptr[j*ldb];
        }
      }
    }
  }

  inline void dcsrmmTT(const int m, const int n, const int k, const double alpha,
    const double *A, const int *rowptr, const int *colind, const int nnz, const int origin,
    const double *B, const int ldb, const double beta, double *C, const int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = &B[0];
        double* Cptr = &C[0];
        for (int v=0; v<k; v++, Bptr+=1, Cptr+=1) {
          Cptr[j*ldc] += alpha * A[z] * Bptr[i*ldb];
        }
      }
    }
  }

  //   template <typename ValueT>
  //   inline
  //   dense_matrix<ValueT>& dgemm(
  //     const trans_t& transA,
  //     const trans_t& transB,
  //     const ValueT alpha,
  //     const csr_matrix<ValueT>& A,
  //     const dense_matrix<ValueT>& B,
  //     const ValueT beta,
  //     dense_matrix<ValueT>& C) {
  //
  //     switch (transB) {
  //       case NoTrans:
  //       switch (transA) {
  //         case NoTrans:
  //         C *= beta;
  //         for (size_type i=0; i<A.nrow(); i++) {
  //           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
  //             size_type j = A.colind(z);
  //             for (index_type jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
  //               C(i + C.rbegin(), jy) += alpha * A.value(z) * B(j + B.rbegin(), jx);
  //             }
  //           }
  //         }
  //         break;
  //         case Trans:
  //         C *= beta;
  //         for (size_type i=0; i<A.nrow(); i++) {
  //           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
  //             size_type j = A.colind(z);
  //             for (index_type jx=B.cbegin(), jy=C.cbegin(); jx<=B.cend(); jx++, jy++) {
  //               C(j + C.rbegin(), jy) += alpha * A.value(z) * B(i + B.rbegin(), jx);
  //             }
  //           }
  //         }
  //         break;
  //       }
  //       break;
  //       case Trans:
  //       switch (transA) {
  //         case NoTrans:
  //         C *= beta;
  //         for (size_type i=0; i<A.nrow(); i++) {
  //           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
  //             size_type j = A.colind(z);
  //             for (index_type jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
  //               C(jy, i + C.cbegin()) += alpha * A.value(z) * B(jx, j + B.cbegin());
  //             }
  //           }
  //         }
  //         break;
  //         case Trans:
  //         C *= beta;
  //         for (size_type i=0; i<A.nrow(); i++) {
  //           for (size_type z=A.rowptr(i); z<A.rowptr(i+1); z++) {
  //             size_type j = A.colind(z);
  //             for (index_type jx=B.rbegin(), jy=C.rbegin(); jx<=B.rend(); jx++, jy++) {
  //               C(jy, j + C.cbegin()) += alpha * A.value(z) * C(jx, i + C.cbegin());
  //             }
  //           }
  //         }
  //         break;
  //       }
  //       break;
  //     }
  //     return C;
  //   }
  //

}
