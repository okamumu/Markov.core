/*
  test
*/

#include <iostream>
#include <iomanip>
#include <random>

#include "marlib.h"

#include "_vector4test.h"
#include "_dense_matrix4test.h"
#include "_csr_matrix4test.h"

// #include "cppblas_csr.hpp"
//
// #include "csc_matrix.h"
// #include "cppblas_csc.hpp"
//
// #include "coo_matrix.h"
// #include "cppblas_coo.hpp"

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

constexpr int SIZE = 20;

void test_copy() {
  double a[SIZE];
  double b[SIZE];
  double aans[SIZE];
  double bans[SIZE];
  int rowptr[SIZE];
  int colind[SIZE];
  double value[SIZE];
  random_vector(SIZE, a);
  random_vector(SIZE, b);

  int m = 3;
  int n = 2;
  int ld = 4;
  vector va(SIZE, a);
  vector vb(SIZE, b);
  vector vaans(SIZE, aans);
  vector vbans(SIZE, bans);
  dense_matrix<> A(3,2,a);
  dense_matrix<> B(3,2,b);
  dense_matrix<std::false_type> Ad(m,n,a,ld);
  dense_matrix<std::false_type> Bd(m,n,b,ld);

  dcopy(1, va);
  for (int i=0; i<SIZE; i++) {
    aans[i] = 1;
  }
  print(va);
  std::cout << check_equal("copy an integer to a vector", SIZE, a, 1, aans, 1) << std::endl;

  dcopy(vb, va);
  for (int i=0; i<SIZE; i++) {
    aans[i] = b[i];
  }
  print(va);
  std::cout << check_equal("copy a vector to a vector", SIZE, a, 1, aans, 1) << std::endl;

  dcopy(0.0, va);
  for (int i=0; i<SIZE; i++) {
    aans[i] = 0.0;
  }
  print(va);
  std::cout << check_equal("copy a double to a vector", SIZE, a, 1, aans, 1) << std::endl;

  dcopy(1, A);
  for (int j=0; j<n; j++) {
    double* ptr = &aans[j*m];
    for (int i=0; i<m; i++) {
      ptr[i] = 1;
    }
  }
  print(va);
  print(vaans);
  std::cout << check_equal("copy a double to a matrix", SIZE, a, 1, aans, 1) << std::endl;

  dcopy(0, va);
  dcopy(0, vaans);
  dcopy(1, Ad);
  for (int j=0; j<n; j++) {
    double* ptr = &aans[j*ld];
    for (int i=0; i<m; i++) {
      ptr[i] = 1;
    }
  }
  print(va);
  print(vaans);
  std::cout << check_equal("copy a double to a matrix (ld)", SIZE, a, 1, aans, 1) << std::endl;

  dcopy(vb, va);
  dcopy(vb, vaans);
  dcopy(0, vb);
  dcopy(0, vbans);
  dcopy(Ad, B);
  for (int j=0; j<n; j++) {
    double* ptr = &aans[j*ld];
    double* bptr = &bans[j*m];
    for (int i=0; i<m; i++) {
      bptr[i] = ptr[i];
    }
  }
  print(vb);
  print(vbans);
  std::cout << check_equal("copy a vector/matrix to a matrix (ld)", SIZE, b, 1, bans, 1) << std::endl;

  {
    random_vector(SIZE, a);
    int nnz = dnnz(A);
    int origin = 1;
    csr_matrix spA(m, n, rowptr, colind, value, nnz, origin);
    dcopy(A, spA);
    print(A);
    print(spA);
    dcopy(0, va);
    dcopy(0, vaans);
    va[0] = nnz;
    vaans[0] = dnnz(spA);
    std::cout << check_equal("nnz for csr", SIZE, a, 1, aans, 1) << std::endl;
  }

  {
    random_vector(SIZE, b);
    dcopy(vb, va);
    dcopy(vb, vaans);
    int nnz = dnnz(A);
    int origin = 1;
    csr_matrix spA(m, n, rowptr, colind, value, nnz, origin);
    dcopy(A, spA);
    dcopy(eye_matrix(), A);
    dcopy(eye_matrix(), spA);
    print(A);
    print(spA);
    dense_matrix<> Aans(3,2,aans);
    dcopy(spA, Aans);
    print(Aans);
    for (int i=0; i<SIZE; i++) {
      va[i] *= vb[i];
    }
    // mapapply(va, [=](double& v) { v *= vb[i]; });
    for (int i=0; i<SIZE; i++) {
      vaans[i] *= vb[i];
    }
    print(va);
    print(vaans);
    std::cout << check_equal("sparse (csr)", SIZE, a, 1, aans, 1) << std::endl;
  }
}

// void test_plus() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(va);
//
//   dplusa(va, 10.0);
//   print(va);
//   for (int i=0; i<SIZE; i++) {
//     c[i] += 10;
//   }
//   print(vc);
//   std::cout << check_equal("plus a", SIZE, c, 1, a, 1) << std::endl;
//
//   dmulti(va, vb);
//   print(va);
//   for (int i=0; i<SIZE; i++) {
//     c[i] *= b[i];
//   }
//   print(vc);
//   std::cout << check_equal("multi", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_csrN() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csr_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMN(A, vb);
//   print(A);
//
//   dmultiMN(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM csr n", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_csrT() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csr_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMT(A, vb);
//   print(A);
//
//   dmultiMT(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM csr t", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_cscN() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csc_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMN(A, vb);
//   print(A);
//
//   dmultiMN(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM csc n", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_cscT() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csc_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMT(A, vb);
//   print(A);
//
//   dmultiMT(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM csc t", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_cooN() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   coo_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMN(A, vb);
//   print(A);
//
//   dmultiMN(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM coo n", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_plus_cooT() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   dense_matrix C(m, n, c);
//   vector va(SIZE, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   dcopy(va, vc);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   coo_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   dmultiMT(A, vb);
//   print(A);
//
//   dmultiMT(spA, vb);
//   print(spA);
//   print(vc);
//   dcopy(spA, C);
//
//   print(va);
//   print(vc);
//   std::cout << check_equal("dmultiM coo t", SIZE, c, 1, a, 1) << std::endl;
// }
//
// void test_mv_csr() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csr_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   double alpha = dist(engine);
//   double beta = dist(engine);
//   double x[SIZE];
//   double y1[SIZE];
//   double y2[SIZE];
//
//   vector vx(SIZE, x);
//   vector vy1(SIZE, y1);
//   vector vy2(SIZE, y2);
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvN(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvN(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv csr n", SIZE, y1, 1, y2, 1) << std::endl;
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvT(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvT(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv csr t", SIZE, y1, 1, y2, 1) << std::endl;
// }
//
// void test_mv_csc() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   csc_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   double alpha = dist(engine);
//   double beta = dist(engine);
//   double x[SIZE];
//   double y1[SIZE];
//   double y2[SIZE];
//
//   vector vx(SIZE, x);
//   vector vy1(SIZE, y1);
//   vector vy2(SIZE, y2);
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvN(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvN(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv csc n", SIZE, y1, 1, y2, 1) << std::endl;
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvT(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvT(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv csc t", SIZE, y1, 1, y2, 1) << std::endl;
// }
//
// void test_mv_coo() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//   coo_matrix spA(m, n, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   double alpha = dist(engine);
//   double beta = dist(engine);
//   double x[SIZE];
//   double y1[SIZE];
//   double y2[SIZE];
//
//   vector vx(SIZE, x);
//   vector vy1(SIZE, y1);
//   vector vy2(SIZE, y2);
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvN(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvN(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv coo n", SIZE, y1, 1, y2, 1) << std::endl;
//
//   dcopy(vb, vx);
//   dcopy(vc, vy1);
//   dgemvT(alpha, A, vx, beta, vy1);
//   print(vy1);
//   dcopy(vc, vy2);
//   dgemvT(alpha, spA, vx, beta, vy2);
//   print(vy2);
//   std::cout << check_equal("dgemv coo t", SIZE, y1, 1, y2, 1) << std::endl;
// }
//
// void test_dgemm_csr1() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   int k = 5;
//
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//
//   dense_matrix A(m, k, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//
//   csr_matrix spA(m, k, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   double alpha = dist(engine);
//   double beta = dist(engine);
//   double x[SIZE];
//   double y1[SIZE];
//   double y2[SIZE];
//
//   {
//     dblas::dcopy(SIZE, b, 1, x, 1);
//     dblas::dcopy(SIZE, c, 1, y1, 1);
//     dblas::dcopy(SIZE, c, 1, y2, 1);
//     dense_matrix vx(k, n, x);
//     dense_matrix vy1(m, n, y1);
//     dense_matrix vy2(m, n, y2);
//     marlib::dgemm('N', 'N', alpha, A, vx, beta, vy1);
//     print(vy1);
//     marlib::dgemm('N', 'N', alpha, spA, vx, beta, vy2);
//     print(vy2);
//     std::cout << check_equal("dgemm csr n n", SIZE, y1, 1, y2, 1) << std::endl;
//   }
//
//   {
//     dblas::dcopy(SIZE, b, 1, x, 1);
//     dblas::dcopy(SIZE, c, 1, y1, 1);
//     dblas::dcopy(SIZE, c, 1, y2, 1);
//     dense_matrix vx(n, k, x);
//     dense_matrix vy1(m, n, y1);
//     dense_matrix vy2(m, n, y2);
//     marlib::dgemm('N', 'T', alpha, A, vx, beta, vy1);
//     print(vy1);
//     marlib::dgemm('N', 'T', alpha, spA, vx, beta, vy2);
//     print(vy2);
//     std::cout << check_equal("dgemm csr n t", SIZE, y1, 1, y2, 1) << std::endl;
//   }
// }
//
// void test_dgemm_csr2() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 2;
//   int k = 5;
//
//   double spa[SIZE];
//   int rowptr[SIZE];
//   int colind[SIZE];
//
//   dense_matrix A(k, m, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   int nz = dnnz(A);
//   int origin = 1;
//
//   csr_matrix spA(k, m, rowptr, colind, spa, nz, origin);
//   dcopy(A, spA);
//   print(spA);
//
//   double alpha = dist(engine);
//   double beta = dist(engine);
//   double x[SIZE];
//   double y1[SIZE];
//   double y2[SIZE];
//
//   {
//     dblas::dcopy(SIZE, b, 1, x, 1);
//     dblas::dcopy(SIZE, c, 1, y1, 1);
//     dblas::dcopy(SIZE, c, 1, y2, 1);
//     dense_matrix vx(k, n, x);
//     dense_matrix vy1(m, n, y1);
//     dense_matrix vy2(m, n, y2);
//     marlib::dgemm('T', 'N', alpha, A, vx, beta, vy1);
//     print(vy1);
//     marlib::dgemm('T', 'N', alpha, spA, vx, beta, vy2);
//     print(vy2);
//     std::cout << check_equal("dgemm csr t n", SIZE, y1, 1, y2, 1) << std::endl;
//   }
//
//   {
//     dblas::dcopy(SIZE, b, 1, x, 1);
//     dblas::dcopy(SIZE, c, 1, y1, 1);
//     dblas::dcopy(SIZE, c, 1, y2, 1);
//     dense_matrix vx(n, k, x);
//     dense_matrix vy1(m, n, y1);
//     dense_matrix vy2(m, n, y2);
//     marlib::dgemm('T', 'T', alpha, A, vx, beta, vy1);
//     print(vy1);
//     marlib::dgemm('T', 'T', alpha, spA, vx, beta, vy2);
//     print(vy2);
//     std::cout << check_equal("dgemm csr t t", SIZE, y1, 1, y2, 1) << std::endl;
//   }
// }
//
// void test_eye() {
//   double a[SIZE];
//   double b[SIZE];
//   double c[SIZE];
//   random_vector(SIZE, a);
//   random_vector(SIZE, b);
//   random_vector(SIZE, c);
//
//   int m = 3;
//   int n = 3;
//   // int ld = 4;
//   dense_matrix A(m, n, a);
//   vector vb(SIZE, b);
//   vector vc(SIZE, c);
//   print(A);
//
//   dfill(A, eye());
//
//   std::cout << "test eye" << std::endl;
//   std::cout << A << std::endl;
//
// }

int main() {
  test_copy();

  // test_plus();
  // test_plus_csrN();
  // test_plus_csrT();
  // test_plus_cscN();
  // test_plus_cscT();
  // test_plus_cooN();
  // test_plus_cooT();
  //
  // test_mv_csr();
  // test_mv_csc();
  // test_mv_coo();
  //
  // test_dgemm_csr1();
  // test_dgemm_csr2();
  //
  // test_eye();
}
