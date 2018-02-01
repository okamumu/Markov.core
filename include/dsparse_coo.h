
namespace dblas {

  inline
  void dense_to_coo(int m, int n, const double *A, int lda,
    double *spA, int *rowind, int *colind, int nnz, int origin) {
    for (int j=0; j<n; j++) {
      const double *Aptr = A + j*lda;
      for (int i=0; i<m; i++, Aptr+=1) {
        if (*Aptr != 0) {
          *spA = *Aptr;
          spA++;
          *rowind = i + origin;
          rowind++;
          *colind = j + origin;
          colind++;
        }
      }
    }
  }

  inline
  void coo_to_dense(int m, int n,
    const double *spA, const int *rowind, const int *colind, int nnz, int origin,
    double *A, int lda) {
    dfill(m, n, A, lda, 0);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      A[i+j*lda] = spA[z];
    }
  }

  // level 2

  inline void dcoomvN(int m, int n, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(m, beta, y, incy);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      y[i*incy] += alpha * A[z] * x[j*incx];
    }
  }

  inline void dcoomvT(int m, int n, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *x, int incx, double beta, double *y, int incy) {
    dscal(n, beta, y, incy);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      y[j*incy] += alpha * A[z] * x[i*incx];
    }
  }

  inline void dcoor(int m, int n, double alpha,
    const double *x, int incx, double *y, int incy,
    double *A, const int *rowind, const int *colind, int nnz, int origin) {
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      A[z] += alpha * x[i*incx] * y[j*incy];
    }
  }

  inline void dcoommNN(int m, int n, int k, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      const double* Bptr = B;
      double* Cptr = C;
      for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
        Cptr[i] += alpha * A[z] * Bptr[j];
      }
    }
  }

  inline void dcoommTN(int m, int n, int k, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      const double* Bptr = B;
      double* Cptr = C;
      for (int v=0; v<n; v++, Bptr+=ldb, Cptr+=ldc) {
        Cptr[j] += alpha * A[z] * Bptr[i];
      }
    }
  }

  inline void dcoommNT(int m, int n, int k, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      const double* Bptr = B;
      double* Cptr = C;
      for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
        Cptr[i] += alpha * A[z] * Bptr[j*ldb];
      }
    }
  }

  inline void dcoommTT(int m, int n, int k, double alpha,
    const double *A, const int *rowind, const int *colind, int nnz, int origin,
    const double *B, int ldb, double beta, double *C, int ldc) {
    dscal(m, n, beta, C, ldc);
    for (int z=0; z<nnz; z++) {
      int i = rowind[z] - origin;
      int j = colind[z] - origin;
      const double* Bptr = B;
      double* Cptr = C;
      for (int v=0; v<n; v++, Bptr+=1, Cptr+=ldc) {
        Cptr[j] += alpha * A[z] * Bptr[i*ldb];
      }
    }
  }
}
