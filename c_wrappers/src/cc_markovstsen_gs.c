///

void f90_markovstsen_gs_dense_(int* n,
  double* Q, int* ldq, double* dQ, int* lddq,
  double* pis, int* incp, double* sstart, int* incss,
  double* s, int* incs, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovstsen_gs_csr_(int* n,
  double* Q, int* rowptr0, int* colind0, int* nnz0,
  double* dQ, int* rowptr1, int* colind1, int* nnz1,
  double* pis, int* incp, double* sstart, int* incss,
  double* s, int* incs, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

void f90_markovstsen_gs_csc_(int* n,
  double* Q, int* colptr0, int* rowind0, int* nnz0,
  double* dQ, int* colptr1, int* rowind1, int* nnz1,
  double* pis, int* incp, double* sstart, int* incss,
  double* s, int* incs, int* maxiter, double* rtol, int* steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*));

///

void cc_markovstsen_gs_dense(int n,
  double* Q, int ldq, double* dQ, int lddq,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_dense_(&n, Q, &ldq, dQ, &lddq, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovstsen_gs_csr(int n,
  double* Q, int* rowptr0, int* colind0, int nnz0,
  double* dQ, int* rowptr1, int* colind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_csr_(&n, Q, rowptr0, colind0, &nnz0,
    dQ, rowptr1, colind1, &nnz1, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}

void cc_markovstsen_gs_csc(int n,
  double* Q, int* colptr0, int* rowind0, int nnz0,
  double* dQ, int* colptr1, int* rowind1, int nnz1,
  double* pis, int incp, double* sstart, int incss,
  double* s, int incs, int maxiter, double rtol, int steps,
  int* iter, double* rerror, int* info, void (*callback)(int*,double*)) {

  f90_markovstsen_gs_csc_(&n, Q, colptr0, rowind0, &nnz0,
    dQ, colptr1, rowind1, &nnz1, pis, &incp,
    sstart, &incss, s, &incs, &maxiter, &rtol, &steps,
    iter, rerror, info, callback);
}
