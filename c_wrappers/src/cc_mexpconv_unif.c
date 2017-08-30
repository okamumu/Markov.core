///

void f90_mexpconv_unif_dense_vec_(char* transQ, char* transH,
  int* n, double* P, int* ldp, double* qv,
  int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* z, int* incz,
  double* H, int* ldh, double* atol);

void f90_mexpconv_unif_csr_vec_(char* transQ, char* transH,
  int* n, double* spP, int* rowptr, int* colind, int* nnz,
  double* qv, int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* z, int* incz,
  double* spH, int* atol);

void f90_mexpconv_unif_csc_vec_(char* transQ, char* transH,
  int* n, double* spP, int* colptr, int* rowind, int* nnz,
  double* qv, int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* z, int* incz,
  double* spH, int* atol);

void f90_mexpconv_unif_coo_vec_(char* transQ, char* transH,
  int* n, double* spP, int* rowind, int* colind, int* nnz,
  double* qv, int* left, int* right, double* poi, double* weight,
  double* x, int* incx, double* y, int* incy, double* z, int* incz,
  double* spH, int* atol);

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
