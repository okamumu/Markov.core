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
  double spa[nz];
  int rowptr[m+1];
  int colind[nz];
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

void test_mm2() {
  double A_[] = {1,2,3,4,5,6,7,8,9,10,11,12};
  dense_matrix A(4,3,A_);
  std::cout << A << std::endl;

  double B_[] = {1,2,3,3,1,4,6,4,2,1,3,4};
  dense_matrix B(4,3, B_);
  std::cout << B << std::endl;

  double C_[9];
  dense_matrix C1(3,3,C_);
  dfill(C1, 0);
  std::cout << C1 << std::endl;

  std::cout << dgemm('T', 'N', 1, A, B, 0, C1) << std::endl;
}

void test_mm3() {
  double A_[] = {1,2,3,4,5,6,7,8,9,10,11,12};
  dense_matrix A(3,4,A_);
  std::cout << A << std::endl;

  double B_[] = {1,2,3,3,1,4,6,4,2,1,3,4};
  dense_matrix B(3,4, B_);
  std::cout << B << std::endl;

  double C_[9];
  dense_matrix C1(3,3,C_);
  dfill(C1, 0);
  std::cout << C1 << std::endl;

  std::cout << dgemm('N', 'T', 1, A, B, 0, C1) << std::endl;
}

void test_mm4() {
  double A_[] = {1,2,3,4,5,6,7,8,9,10,11,12};
  dense_matrix A(4,3,A_);
  std::cout << A << std::endl;

  double B_[] = {1,2,3,3,1,4,6,4,2,1,3,4};
  dense_matrix B(3,4, B_);
  std::cout << B << std::endl;

  double C_[9];
  dense_matrix C1(3,3,C_);
  dfill(C1, 0);
  std::cout << C1 << std::endl;

  std::cout << dgemm('T', 'T', 1, A, B, 0, C1) << std::endl;
}


// void test_ger1() {
//   vector<double,int> x = {1,2,3};
//   vector<double,int> y = {1,2,3,4};
//
//   dense_matrix<double,int> C1(3,4);
//   C1 = 0;
//   // std::cout << C1.ger(NoTrans, 10, x, y) << std::endl;
//   std::cout << dger<double,int>(NoTrans, 10, x, y, C1) << std::endl;
// }
//
// void test_ger2() {
//   vector<double,int> x = {1,2,3};
//   vector<double,int> y = {1,2,3,4};
//
//   dense_matrix<double,int> C1(4,3);
//   C1 = 0;
//   // std::cout << C1.ger(Trans, 10, x, y) << std::endl;
//   std::cout << dger<double,int>(Trans, 10, x, y, C1) << std::endl;
// }
//
// void test_gesv1() {
//   std::cout << "gesv1" << std::endl;
//   dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
//   std::cout << A << std::endl;
//
//   dense_matrix<double,int> B(3,3, {1,0,0,0,1,0,0,0,1});
//   std::cout << B << std::endl;
//
//   dense_matrix<double,int> C1(3,3);
//   C1 = 0;
//   std::cout << C1 << std::endl;
//
//   // std::cout << C1.gesv(A, B) << std::endl;
//   C1 = B;
//   std::cout << dgesv(A, C1) << std::endl;
//   std::cout << A << std::endl;
//   std::cout << B << std::endl;
//   std::cout << C1 << std::endl;
// }
//
// void test_gesv2() {
//   std::cout << "gesv2" << std::endl;
//   dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
//   std::cout << A << std::endl;
//
//   dense_matrix<double,int> B(3,2, {3,4,5,3,2,1});
//   std::cout << B << std::endl;
//
//   dense_matrix<double,int> C1(3,2);
//   C1 = 0;
//   std::cout << C1 << std::endl;
//
//   // std::cout << C1.gesv(A, B) << std::endl;
//   C1 = B;
//   std::cout << dgesv(A, C1) << std::endl;
//   std::cout << A << std::endl;
//   std::cout << B << std::endl;
//   std::cout << C1 << std::endl;
// }
//
// void test_gesv3() {
//   std::cout << "gesv3" << std::endl;
//   dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
//   std::cout << A << std::endl;
//
//   dense_matrix<double,int> B = dense_matrix<double,int>::eye(3);
//   // B.eye();
//   std::cout << B << std::endl;
//
//   // std::cout << B.gesv(A) << std::endl;
// }
//
// void test_array0() {
//   std::cout << "array 0" << std::endl;
//   vector<double,int> x = {-1,2,3};
//   array<vector<double,int>*> y(3);
//   std::cout << "set1" << std::endl;
//   y.ptr(0) = new vector<double,int>(x.clone());
//   y.ptr(1) = new vector<double,int>(x.clone());
//   y.ptr(2) = new vector<double,int>(x.clone());
//   std::cout << y << std::endl;
//   y[0](1) = 10;
//   std::cout << y << std::endl;
//   array<vector<double,int>*> z = y.subarray(1);
//   std::cout << z << std::endl;
// }
//
// void test_array1() {
//   std::cout << "array 1" << std::endl;
//   vector<double*,int> x(3);
//   for (int i=x.begin(); i<=x.end(); i++) {
//     x.ptr(i) = new double(i*100);
//   }
//   std::cout << x << std::endl;
//
//   vector<double,int> y(3);
//   y = x;
//   std::cout << y << std::endl;
// }
//
// void test_tr1() {
//   std::cout << "tr 1" << std::endl;
//   dense_matrix<double,int> ma(3,4,{1,2,3,4,5,6,7,8,9,10,11,12});
//   std::cout << ma << std::endl;
//   std::cout << ma.tr() << std::endl;
// }
//
// void test_diag1() {
//   std::cout << "diag1" << std::endl;
//   dense_matrix<double,int> A(3,3, {-1,2,3,4,-5,6,7,8,-9});
//   std::cout << A << std::endl;
//
//   vector<double*,int> d(3);
//   double v = 0;
//
//   d.ptr(1) = &v;
//   d.ptr(2) = &v;
//   d.ptr(3) = &v;
//
//   A.diag(d, 1, d.begin());
//
//   std::cout << d << std::endl;
// }
//
// void test0() {
//   // double *a = new double [10];
//   // for (int i=0; i<10; i++) {
//   //   a[i] = i;
//   // }
//   //
//   // std::cout << a[0] << std::endl;
//   //
//   // double& ref_a = a[5];
//   // double* ptr_a = &(a[2]);
//   //
//   // std::cout << ref_a << std::endl;
//   // std::cout << *ptr_a << std::endl;
//   //
//   // delete [] a;
//   // a = NULL;
//   //
//   // std::cout << ref_a << std::endl;
//   // std::cout << *ptr_a << std::endl;
//
//   array<double> a(3);
//   a[0] = 10;
//   a[1] = 1;
//   a[2] = 1000;
//   std::cout << a << std::endl;
//
//   array<double> b = a;
//   std::cout << b << std::endl;
//
//   b[1] = 10;
//   std::cout << a << std::endl;
//   std::cout << b << std::endl;
//
//   range<int> r(1,2);
//
//   // vector<double,int> v(range<int>(-1,10));
//   vector<double,int> v(12);
//   for (int i=v.begin(), x=10; i<=v.end(); i++, x++) {
//     v(i) = x;
//     std::cout << i << " " << v(i) << std::endl;
//   }
//   std::cout << v.value() << std::endl;
//
//   const vector<double,int> y = v(range<int>(2,10));
//   for (int i=y.begin(); i<=y.end(); i++) {
//     std::cout << i << " " << y(i) << std::endl;
//   }
//   std::cout << y.value() << std::endl;
//
//   std::cout << y(range<int>(3,5)) << std::endl;
//
//   v.set_range(range<int>(1,12));
//   std::cout << "v" << std::endl;
//   for (int i=v.begin(); i<=v.end(); i++) {
//     std::cout << i << " " << v(i) << std::endl;
//   }
//
//   vector<double,int> w = y.clone();
//   std::cout << "w" << std::endl;
//   for (int i=w.begin(); i<=w.end(); i++) {
//     std::cout << i << " " << w(i) << std::endl;
//   }
//
//   v = 100;
//   for (int i=v.begin(); i<=v.end(); i++) {
//     std::cout << i << " " << v(i) << std::endl;
//   }
//
//   w += y;
//   for (int i=w.begin(); i<=w.end(); i++) {
//     std::cout << i << " " << w(i) << std::endl;
//   }
//
//   vector<double,int> z = y + w;
//   for (int i=z.begin(); i<=z.end(); i++) {
//     std::cout << i << " " << z(i) << std::endl;
//   }
//
//   std::cout << z << std::endl;
//
//   std::cout << iamax(z) << std::endl;
//
//   std::cout << "matrix" << std::endl;
//
//   dense_matrix<double,int> m1(10,12);
//
//   m1 = 1;
//   std::cout << m1 << std::endl;
//
//   for (int j=m1.cbegin(); j<=m1.cend(); j++) {
//     for (int i=m1.rbegin(); i<=m1.rend(); i++) {
//       m1(i,j) = 13 * i + j;
//     }
//   }
//
//   const dense_matrix<double,int> m2 = m1;
//   std::cout << m2 << std::endl;
//
//   dense_matrix<double,int> m3 = m1.clone();
//   std::cout << m3 << std::endl;
//
//   std::cout << m3(range<int>(1,5),range<int>(1,3)) << std::endl;
//
//   std::cout << m1 + m3 << std::endl;
//   std::cout << m1 - m3 << std::endl;
//   std::cout << m1 * m3 << std::endl;
//   std::cout << m1 / m3 << std::endl;
//
//   std::cout << m1 << std::endl;
//
//   std::cout << "matrix" << std::endl;
//   vector<double,int> xv = {1,2,3,4,5,6,7,8,9,10,11,12};
//   std::cout << xv << std::endl;
//
//   vector<double,int> yv(10);
//   yv = 0;
//   std::cout << yv << std::endl;
//
//   // yv.gemv(NoTrans, 1, m1, xv, 0);
//   dgemv<double,int>(NoTrans, 1, m1, xv, 0, yv);
//   std::cout << yv << std::endl;
//
// }

int main() {
  test_mv_csr();
  test_mm2();
  test_mm3();
  test_mm4();
}
