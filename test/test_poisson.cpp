/*
  test
*/

#include <iostream>
#include <iomanip>
#include <random>

#include "marlib.h"

#include "vector.h"

using namespace marlib;

std::random_device seed_gen;
std::mt19937 engine(seed_gen());
std::uniform_real_distribution<> dist(0, 10000.0);

std::string check_equal(const std::string& title, int n, const double *x, int incx, const double *y, int incy) {
  std::cout << title << ": ";
  for (int i=0; i<n; i++, x+=incx, y+=incy) {
    if (std::abs(*x - *y) > 1.0e-12) {
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

void test_poi1() {
  std::cout << "poi1" << std::endl;

  double lambda = dist(engine);
  int right = poisson<double>::rightbound(lambda, 1.0e-8);
  double prob[right+1];
  vector vprob(right+1, prob);
  poisson<double> poi(right+1, prob);

  poi.set(lambda, 0, right);
  print(vprob);

  double sum = 0;
  double one = 1;
  for (int i=poi.left(); i<=poi.right(); i++) {
    sum += poi(i);
  }
  sum /= poi.weight();
  std::cout << check_equal("test poi", 1, &sum, 1, &one, 1) << std::endl;
}

int main() {
  test_poi1();
}
