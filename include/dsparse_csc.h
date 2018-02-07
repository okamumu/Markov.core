namespace dblas {

  inline
  void dense_to_csc(int m, int n, const double *A, int lda,
    double *spA, int *colptr, int *rowind, int, int origin) {
    int z = origin;
    for (int j=0; j<n; j++) {
      *colptr = z;
      colptr++;
      const double *Aptr = A + j*lda;
      for (int i=0; i<m; i++, Aptr+=1) {
        if (!is_zero(*Aptr)) {
          *spA = *Aptr;
          spA++;
          *rowind = i + origin;
          rowind++;
          z++;
        }
      }
    }
    *colptr = z;
  }

  inline
  void csc_to_dense(int m, int n,
    const double *spA, const int *colptr, const int *rowind, int, int origin,
    double *A, int lda) {
    dfill(m, n, A, lda, 0);
    for (int j=0; j<n; j++, A+=lda) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        A[i] = spA[z];
      }
    }
  }

  // level 2

  inline void dcscmvN(int m, int n, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(m, beta, y, incy);
    for (int j=0; j<n; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        y[i*incy] += alpha * A[z] * x[j*incx];
      }
    }
  }

  inline void dcscmvT(int, int n, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(n, beta, y, incy);
    for (int j=0; j<n; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        y[j*incy] += alpha * A[z] * x[i*incx];
      }
    }
  }

  inline void dcscr(int, int n, double alpha,
    const double *x, int incx, const double *y, int incy,
    double *A, const int *colptr, const int *rowind, int, int origin) {
    for (int j=0; j<n; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        A[z] += alpha * x[i*incx] * y[j*incy];
      }
    }
  }

  inline void dcscmmNN(int m, int n, int k, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int j=0; j<k; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j];
        }
      }
    }
  }

  inline void dcscmmTN(int m, int n, int, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int j=0; j<m; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i];
        }
      }
    }
  }

  inline void dcscmmNT(int m, int n, int k, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int j=0; j<k; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[i] += alpha * A[z] * Bptr[j*ldb];
        }
      }
    }
  }

  inline void dcscmmTT(int m, int n, int, double alpha,
    const double *A, const int *colptr, const int *rowind, int, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int j=0; j<m; j++) {
      for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
        int i = rowind[z] - origin;
        const double* Bptr = B;
        double* Cptr = C;
        for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
          Cptr[j] += alpha * A[z] * Bptr[i*ldb];
        }
      }
    }
  }

}
