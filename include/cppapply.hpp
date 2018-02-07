/*
  blas.hpp
*/

namespace marlib {

  template <typename VectorT, typename Func>
  VectorT& map_vector(VectorT& x, Func f) {
    for (int i=0; i<size(x); i++) {
      f(i, x[i]);
    }
    return x;
  }

  template <typename MatrixT, typename Func>
  MatrixT& map_matrix(MatrixT& x, Func f) {
    double *xptr = &x[0];
    for (int j=0; j<ncol(x); j++, xptr+=ld(x)) {
      for (int i=0; i<nrow(x); i++) {
        f(i, j, xptr[i]);
      }
    }
    return x;
  }

  ////

  template <typename VectorT>
  VectorT& dplusa(VectorT& x, double a) {
    return map_vector(x, [=](int, double& value){ value += a; });
  }

  template <typename VectorT>
  VectorT& dmulti(VectorT& x, const VectorT& y) {
    return map_vector(x, [=](int i, double& value){ value *= y[i]; });
  }

  template <typename VectorT>
  VectorT& ddiv(VectorT& x, const VectorT& y) {
    return map_vector(x, [=](int i, double& value){ value /= y[i]; });
  }

  template <typename MatrixT, typename VectorT>
  MatrixT& dmultiMN(MatrixT& x, const VectorT& y) {
    return map_matrix(x, [=](int i, int, double& value){ value *= y[i]; });
  }

  template <typename MatrixT, typename VectorT>
  MatrixT& dmultiMT(MatrixT& x, const VectorT& y) {
    return map_matrix(x, [=](int, int j, double& value){ value *= y[j]; });
  }

  template <typename MatrixT, typename VectorT>
  MatrixT& ddivMN(MatrixT& x, const VectorT& y) {
    return map_matrix(x, [=](int i, int, double& value){ value /= y[i]; });
  }

  template <typename MatrixT, typename VectorT>
  MatrixT& ddivMT(MatrixT& x, const VectorT& y) {
    return map_matrix(x, [=](int, int j, double& value){ value /= y[j]; });
  }
}
