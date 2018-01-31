
namespace dblas {

  inline
  void dense_to_csr(int m, int n, const double *A, int lda,
    double *spA, int *rowptr, int *colind, int nnz, int origin) {
    int z = origin;
    for (int i=0; i<m; i++) {
      *rowptr = z;
      rowptr++;
      const double *Aptr = A + i;
      for (int j=0; j<n; j++, Aptr+=lda) {
        if (*Aptr != 0) {
          *spA = *Aptr;
          spA++;
          *colind = j + origin;
          colind++;
          z++;
        }
      }
    }
    *rowptr = z;
  }

  inline
  void csr_to_dense(int m, int n,
    const double *spA, const int *rowptr, const int *colind, int nnz, int origin,
    double *A, int lda) {
    double* Aptr = A;
    for (int j=0; j<n; j++, Aptr+=lda) {
      for (int i=0; i<m; i++) {
        Aptr[i] = 0;
      }
    }
    for (int i=0; i<m; i++, A+=1) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        A[j*lda] = spA[z];
      }
    }
  }

  // level 2

  inline void dcsrmvN(int m, int n, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(m, beta, y, incy);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        y[i*incy] += alpha * A[z] * x[j*incx];
      }
    }
  }

  inline void dcsrmvT(int m, int n, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(n, beta, y, incy);
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        y[j*incy] += alpha * A[z] * x[i*incx];
      }
    }
  }

  inline void dcsrr(int m, int n, double alpha,
    const double *x, int incx, const double *y, int incy,
    double *A, const int *rowptr, const int *colind, int nnz, int origin) {
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        A[z] += alpha * x[i*incx] * y[j*incy];
      }
    }
  }

  inline void dcsrmmNN(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j];
        }
      }
    }
  }

  inline void dcsrmmTN(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i];
        }
      }
    }
  }

  inline void dcsrmmNT(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<m; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j*ldb];
        }
      }
    }
  }

  inline void dcsrmmTT(int m, int n, int k, double alpha,
    const double *A, const int *rowptr, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    double* Cptr = C;
    for (int i=0; i<n; i++, Cptr+=ldc) {
      dscal(m, beta, Cptr, 1);
    }
    for (int i=0; i<k; i++) {
      for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
        int j = colind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i*ldb];
        }
      }
    }
  }

}
