/*
  test
*/

#include <iostream>
#include <iomanip>
#include <random>

#include "dblas.h"
#include "dsparse_csr.h"
#include "dsparse_csc.h"
#include "dsparse_coo.h"

using namespace dblas;

std::random_device seed_gen;
std::mt19937 engine(seed_gen());
std::uniform_real_distribution<> dist(-1.0, 1.0);

std::string check_equal(const std::string& title, int n, const double *x, int incx, const double *y, int incy) {
  std::cout << title << ": ";
  for (int i=0; i<n; i++, x+=incx, y+=incy) {
    if (std::abs(*x - *y) > 1.0e-10) {
      std::cout << "false: " << i << " " << (*x - *y) << std::endl;
      return "false";
    }
  }
  return "true";
}

void print_vector(int n, const double *x, int incx) {
#ifdef DEBUG
  std::cout << "   ";
  for (int i=0; i<n; i++, x+=incx) {
    std::cout << *x << " ";
  }
  std::cout << std::endl;
#endif
}

void print_matrix(int m, int n, const double *x, int ld) {
#ifdef DEBUG
  for (int i=0; i<m; i++) {
    std::cout << "   ";
    for (int j=0; j<n; j++) {
      std::cout << std::setw(10) << x[i+j*ld] << " ";
    }
    std::cout << std::endl;
  }
#endif
}

void print_csr_matrix(int m, int n, const double *A, const int *rowptr, const int *colind, int nnz, int origin) {
#ifdef DEBUG
  for (int i=0; i<m; i++) {
    for (int z=rowptr[i]-origin; z<rowptr[i+1]-origin; z++) {
      int j = colind[z] - origin;
      std::cout << "   ";
      std::cout << "(" << i + origin << "," << j + origin << ")=" << A[z] << std::endl;
    }
  }
#endif
}

void print_csc_matrix(int m, int n, const double *A, const int *colptr, const int *rowind, int nnz, int origin) {
#ifdef DEBUG
  for (int j=0; j<n; j++) {
    for (int z=colptr[j]-origin; z<colptr[j+1]-origin; z++) {
      int i = rowind[z] - origin;
      std::cout << "   ";
      std::cout << "(" << i + origin << "," << j + origin << ")=" << A[z] << std::endl;
    }
  }
#endif
}

void print_coo_matrix(int m, int n, const double *A, const int *colptr, const int *rowind, int nnz, int origin) {
#ifdef DEBUG
    for (int z=0; z<nnz; z++) {
      int i = colptr[z] - origin;
      int j = rowind[z] - origin;
      std::cout << "   ";
      std::cout << "(" << i + origin << "," << j + origin << ")=" << A[z] << std::endl;
    }
#endif
}

void random_vector(int n, double* x) {
  for (int i=0; i<n; i++, x++) {
    *x = dist(engine);
    if (*x < 0) {
      *x = 0;
    }
  }
}

//////////////////////////

void test_dcopy() {
  double a[] = {1,2,4,5,3};
  double b[] = {0,0,0,0,0};
  double c[] = {20, 40};
  double d[] = {20, 2, 40, 5, 3};

  dcopy(5, a, 1, b, 1);
  print_vector(5, a, 1);
  print_vector(5, b, 1);

  std::cout << check_equal("dcopy 1", 5, a, 1, b, 1) << std::endl;

  dcopy(2, c, 1, a, 2);
  print_vector(2, c, 1);
  print_vector(5, a, 1);

  std::cout << check_equal("dcopy 2", 5, a, 1, d, 1) << std::endl;
}

void test_dgemv_csr() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int ld = 4;
  print_matrix(m, n, a, ld);

  int nz = dnnz(m, n, a, ld);
  int origin = 1;
  double spa[nz];
  int rowptr[m+1];
  int colind[nz];
  dense_to_csr(m, n, a, ld, spa, rowptr, colind, nz, origin);
  print_csr_matrix(m, n, spa, rowptr, colind, nz, origin);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csr n", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csr t", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csr n inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csr t inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, a, 1, x, 1);
  csr_to_dense(m, n, spa, rowptr, colind, nz, origin, x, ld);
  print_vector(size, a, 1);
  print_vector(size, x, 1);
  std::cout << check_equal("csr to dense", size, a, 1, x, 1) << std::endl;
}

void test_dgemv_csc() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int ld = 4;
  print_matrix(m, n, a, ld);

  int nz = dnnz(m, n, a, ld);
  int origin = 1;
  double spa[nz];
  int rowptr[n+1];
  int colind[nz];
  dense_to_csc(m, n, a, ld, spa, rowptr, colind, nz, origin);
  print_csc_matrix(m, n, spa, rowptr, colind, nz, origin);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csc n", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csc t", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csc n inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv csc t inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, a, 1, x, 1);
  csc_to_dense(m, n, spa, rowptr, colind, nz, origin, x, ld);
  print_vector(size, a, 1);
  print_vector(size, x, 1);
  std::cout << check_equal("csc to dense", size, a, 1, x, 1) << std::endl;
}

void test_dgemv_coo() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int ld = 4;
  print_matrix(m, n, a, ld);

  int nz = dnnz(m, n, a, ld);
  int origin = 1;
  double spa[nz];
  int rowptr[nz];
  int colind[nz];
  dense_to_coo(m, n, a, ld, spa, rowptr, colind, nz, origin);
  print_coo_matrix(m, n, spa, rowptr, colind, nz, origin);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoomvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv coo n", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 1, beta, y1, 1);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoomvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 1, beta, y2, 1);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv coo t", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('N', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoomvN(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv coo n inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemv('T', m, n, alpha, a, ld, x, 2, beta, y1, 2);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoomvT(m, n, alpha, spa, rowptr, colind, nz, origin, x, 2, beta, y2, 2);
  print_vector(size, y2, 1);
  std::cout << check_equal("dgemv coo t inc=2", size, y1, 1, y2, 1) << std::endl;

  dcopy(size, a, 1, x, 1);
  coo_to_dense(m, n, spa, rowptr, colind, nz, origin, x, ld);
  print_vector(size, a, 1);
  print_vector(size, x, 1);
  std::cout << check_equal("coo to dense", size, a, 1, x, 1) << std::endl;
}

void test_dgemm_csr() {
  int size = 25;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int k = 5;
  int ld = 5;
  int ldb = 5;
  int ldc = 5;

  int nz;
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_csr(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_csr_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmmNN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csr n n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_csr(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_csr_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmmTN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csr t n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_csr(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_csr_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmmNT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csr n t", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_csr(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_csr_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcsrmmTT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csr t t", size, y1, 1, y2, 1) << std::endl;
}

void test_dgemm_csc() {
  int size = 25;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int k = 5;
  int ld = 5;
  int ldb = 5;
  int ldc = 5;

  int nz;
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_csc(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_csc_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmmNN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csc n n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_csc(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_csc_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmmTN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csc t n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_csc(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_csc_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmmNT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csc n t", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_csc(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_csc_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcscmmTT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm csc t t", size, y1, 1, y2, 1) << std::endl;
}

void test_dgemm_coo() {
  int size = 25;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  int k = 5;
  int ld = 5;
  int ldb = 5;
  int ldc = 5;

  int nz;
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_coo(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_coo_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoommNN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm coo n n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_coo(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_coo_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'N', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoommTN(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm coo t n", size, y1, 1, y2, 1) << std::endl;

  print_matrix(m, k, a, ld);
  nz = dnnz(m, k, a, ld);
  dense_to_coo(m, k, a, ld, spa, rowptr, colind, nz, origin);
  print_coo_matrix(m, k, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('N', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoommNT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm coo n t", size, y1, 1, y2, 1) << std::endl;

  print_matrix(k, m, a, ld);
  nz = dnnz(k, m, a, ld);
  dense_to_coo(k, m, a, ld, spa, rowptr, colind, nz, origin);
  print_coo_matrix(k, m, spa, rowptr, colind, nz, origin);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y1, 1);
  dgemm('T', 'T', m, n, k, alpha, a, ld, x, ldb, beta, y1, ldc);
  print_vector(size, y1, 1);
  dcopy(size, b, 1, x, 1);
  dcopy(size, c, 1, y2, 1);
  dcoommTT(m, n, k, alpha, spa, rowptr, colind, nz, origin, x, ldb, beta, y2, ldc);
  print_vector(size, y2, 1);
  print_matrix(ldc, size/ldc, c, ldc);
  print_matrix(ldc, size/ldc, y1, ldc);
  print_matrix(ldc, size/ldc, y2, ldc);
  std::cout << check_equal("dgemm coo t t", size, y1, 1, y2, 1) << std::endl;
}

int main() {
  test_dcopy();
  test_dgemv_csr();
  test_dgemv_csc();
  test_dgemv_coo();

  test_dgemm_csr();
  test_dgemm_csc();
  test_dgemm_coo();
}
