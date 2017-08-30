///

void f90_markovstsen_dgesv_dense_(int* n, double* Q, int* ldq, double* dQ, int* lddq,
  double* pis, int* incp, double* s, int* incs, int* info);

void f90_markovstsen_dgesv_csr_(int* n,
  double* spQ, int* rowptr0, int* colind0, int* nnz0,
  double* dQ, int* rowptr1, int* colind1, int* nnz1,
  double* pis, int* incp, double* s, int* incs, int* info);

void f90_markovstsen_dgesv_csc_(int* n,
  double* spQ, int* colptr0, int* rowind0, int* nnz0,
  double* dQ, int* colptr1, int* rowind1, int* nnz1,
  double* pis, int* incp, double* s, int* incs, int* info);

void f90_markovstsen_dgesv_coo_(int* n,
  double* spQ, int* rowind0, int* colind0, int* nnz0,
  double* dQ, int* rowind1, int* colind1, int* nnz1,
  double* pis, int* incp, double* s, int* incs, int* info);

///

void cc_markovstsen_dgesv_dense(int n, double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_dense_(&n, Q, &ldq, dQ, &lddq, pis, &incp,
    s, &incs, info);
}

void cc_markovstsen_dgesv_csr(int n,
  double* spQ, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_csr_(&n, spQ, rowptr0, colind0, &nnz0,
    dQ, rowptr1, colind1, &nnz1, pis, &incp, s, &incs, info);
}

void cc_markovstsen_dgesv_csc(int n,
  double* spQ, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_csc_(&n, spQ, colptr0, rowind0, &nnz0,
    dQ, colptr1, rowind1, &nnz1, pis, &incp, s, &incs, info);
}

void cc_markovstsen_dgesv_coo(int n,
  double* spQ, int* rowind0, int* colind0, int nnz0,
  double* dQ, int* rowind1, int* colind1, int nnz1,
  double* pis, int incp, double* s, int incs, int* info) {

  f90_markovstsen_dgesv_coo_(&n, spQ, rowind0, colind0, &nnz0,
    dQ, rowind1, colind1, &nnz1, pis, &incp, s, &incs, info);
}
