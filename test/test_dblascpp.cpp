/*
  test
*/

#include <iostream>
#include <iomanip>
#include <random>

#include "marlib.h"

#include "vector.h"
#include "dense_matrix.h"

#include "csr_matrix.h"
#include "cppblas_csr.hpp"

#include "csc_matrix.h"
#include "cppblas_csc.hpp"

#include "coo_matrix.h"
#include "cppblas_coo.hpp"

using namespace marlib;

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

void random_vector(int n, double* x) {
  for (int i=0; i<n; i++, x++) {
    *x = dist(engine);
    if (*x < 0) {
      *x = 0;
    }
  }
}

template <typename T>
void print(const T& x) {
#ifdef DEBUG
  std::cout << x << std::endl;
#endif
}

void test_mv_csr() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  // int ld = 4;
  dense_matrix A(m, n, a);
  vector vb(size, b);
  vector vc(size, c);
  print(A);

  int nz = dnnz(A);
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];
  csr_matrix spA(m, n, rowptr, colind, spa, nz, origin);
  dcopy(A, spA);
  print(spA);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  vector vx(size, x);
  vector vy1(size, y1);
  vector vy2(size, y2);

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvN(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvN(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv csr n", size, y1, 1, y2, 1) << std::endl;

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvT(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvT(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv csr t", size, y1, 1, y2, 1) << std::endl;
}

void test_mv_csc() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  // int ld = 4;
  dense_matrix A(m, n, a);
  vector vb(size, b);
  vector vc(size, c);
  print(A);

  int nz = dnnz(A);
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];
  csc_matrix spA(m, n, rowptr, colind, spa, nz, origin);
  dcopy(A, spA);
  print(spA);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  vector vx(size, x);
  vector vy1(size, y1);
  vector vy2(size, y2);

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvN(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvN(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv csc n", size, y1, 1, y2, 1) << std::endl;

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvT(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvT(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv csc t", size, y1, 1, y2, 1) << std::endl;
}

void test_mv_coo() {
  int size = 20;
  double a[size];
  double b[size];
  double c[size];
  random_vector(size, a);
  random_vector(size, b);
  random_vector(size, c);

  int m = 3;
  int n = 2;
  // int ld = 4;
  dense_matrix A(m, n, a);
  vector vb(size, b);
  vector vc(size, c);
  print(A);

  int nz = dnnz(A);
  int origin = 1;
  double spa[size];
  int rowptr[size];
  int colind[size];
  coo_matrix spA(m, n, rowptr, colind, spa, nz, origin);
  dcopy(A, spA);
  print(spA);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  vector vx(size, x);
  vector vy1(size, y1);
  vector vy2(size, y2);

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvN(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvN(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv coo n", size, y1, 1, y2, 1) << std::endl;

  dcopy(vb, vx);
  dcopy(vc, vy1);
  dgemvT(alpha, A, vx, beta, vy1);
  print(vy1);
  dcopy(vc, vy2);
  dgemvT(alpha, spA, vx, beta, vy2);
  print(vy2);
  std::cout << check_equal("dgemv coo t", size, y1, 1, y2, 1) << std::endl;
}

void test_dgemm_csr1() {
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

  double spa[size];
  int rowptr[size];
  int colind[size];

  dense_matrix A(m, k, a);
  vector vb(size, b);
  vector vc(size, c);
  print(A);

  int nz = dnnz(A);
  int origin = 1;

  csr_matrix spA(m, k, rowptr, colind, spa, nz, origin);
  dcopy(A, spA);
  print(spA);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  {
    dblas::dcopy(size, b, 1, x, 1);
    dblas::dcopy(size, c, 1, y1, 1);
    dblas::dcopy(size, c, 1, y2, 1);
    dense_matrix vx(k, n, x);
    dense_matrix vy1(m, n, y1);
    dense_matrix vy2(m, n, y2);
    marlib::dgemm('N', 'N', alpha, A, vx, beta, vy1);
    print(vy1);
    marlib::dgemm('N', 'N', alpha, spA, vx, beta, vy2);
    print(vy2);
    std::cout << check_equal("dgemm csr n n", size, y1, 1, y2, 1) << std::endl;
  }

  {
    dblas::dcopy(size, b, 1, x, 1);
    dblas::dcopy(size, c, 1, y1, 1);
    dblas::dcopy(size, c, 1, y2, 1);
    dense_matrix vx(n, k, x);
    dense_matrix vy1(m, n, y1);
    dense_matrix vy2(m, n, y2);
    marlib::dgemm('N', 'T', alpha, A, vx, beta, vy1);
    print(vy1);
    marlib::dgemm('N', 'T', alpha, spA, vx, beta, vy2);
    print(vy2);
    std::cout << check_equal("dgemm csr n t", size, y1, 1, y2, 1) << std::endl;
  }
}

void test_dgemm_csr2() {
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

  double spa[size];
  int rowptr[size];
  int colind[size];

  dense_matrix A(k, m, a);
  vector vb(size, b);
  vector vc(size, c);
  print(A);

  int nz = dnnz(A);
  int origin = 1;

  csr_matrix spA(k, m, rowptr, colind, spa, nz, origin);
  dcopy(A, spA);
  print(spA);

  double alpha = dist(engine);
  double beta = dist(engine);
  double x[size];
  double y1[size];
  double y2[size];

  {
    dblas::dcopy(size, b, 1, x, 1);
    dblas::dcopy(size, c, 1, y1, 1);
    dblas::dcopy(size, c, 1, y2, 1);
    dense_matrix vx(k, n, x);
    dense_matrix vy1(m, n, y1);
    dense_matrix vy2(m, n, y2);
    marlib::dgemm('T', 'N', alpha, A, vx, beta, vy1);
    print(vy1);
    marlib::dgemm('T', 'N', alpha, spA, vx, beta, vy2);
    print(vy2);
    std::cout << check_equal("dgemm csr t n", size, y1, 1, y2, 1) << std::endl;
  }

  {
    dblas::dcopy(size, b, 1, x, 1);
    dblas::dcopy(size, c, 1, y1, 1);
    dblas::dcopy(size, c, 1, y2, 1);
    dense_matrix vx(n, k, x);
    dense_matrix vy1(m, n, y1);
    dense_matrix vy2(m, n, y2);
    marlib::dgemm('T', 'T', alpha, A, vx, beta, vy1);
    print(vy1);
    marlib::dgemm('T', 'T', alpha, spA, vx, beta, vy2);
    print(vy2);
    std::cout << check_equal("dgemm csr t t", size, y1, 1, y2, 1) << std::endl;
  }
}

int main() {
  test_mv_csr();
  test_mv_csc();
  test_mv_coo();

  test_dgemm_csr1();
  test_dgemm_csr2();
}
