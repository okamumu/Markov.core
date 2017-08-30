///

void f90_unif_dense_(int* n, double* Q, int* ldq, double* P, int* ldp, double* qv, double* ufact);

void f90_unif_csr_(int* n, double* spQ, int* rowptr, int* colind, int* nnz,
  double* spP, double* qv, double* ufact);

void f90_unif_csc_(int* n, double* spQ, int* colptr, int* rowind, int* nnz,
  double* spP, double* qv, double* ufact);

void f90_unif_coo_(int* n, double* spQ, int* rowind, int* colind, int* nnz,
  double* spP, double* qv, double* ufact);

//

void cc_unif_dense(int n, double* Q, int ldq, double* P, int ldp, double* qv, double ufact) {
  f90_unif_dense_(&n, Q, &ldq, P, &ldp, qv, &ufact);
}

void cc_unif_csr(int n, double* spQ, int* rowptr, int* colind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_csr_(&n, spQ, rowptr, colind, &nnz, spP, qv, &ufact);
}

void cc_unif_csc(int n, double* spQ, int* colptr, int* rowind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_csc_(&n, spQ, colptr, rowind, &nnz, spP, qv, &ufact);
}

void cc_unif_coo(int n, double* spQ, int* rowind, int* colind, int nnz,
  double* spP, double* qv, double ufact) {

  f90_unif_coo_(&n, spQ, rowind, colind, &nnz, spP, qv, &ufact);
}
