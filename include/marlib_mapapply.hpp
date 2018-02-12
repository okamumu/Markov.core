/*
  blas.hpp
*/

namespace marlib {

  template <typename T, typename Func>
  T& mapapply(T& x, Func f) {
    return mapapply_impl(x, f, typename get_category<T>::type());
  }

  template <typename T, typename Func>
  T& mapapply_impl(T& x, Func f, double_vector_tag) {
    using traits = vector_traits<T>;
    auto xptr = traits::begin(x);
    for (int i=0; i<traits::size(x); i++, xptr+=traits::inc(x)) {
      f(*xptr);
    }
    return x;
  }

  template <typename T, typename Func>
  T& mapapply_impl(T& x, Func f, double_dense_matrix_tag) {
    using traits = dense_matrix_traits<T>;
    auto xptr = traits::begin(x);
    for (int j=0; j<traits::ncol(x); j++, xptr+=traits::ld(x)) {
      for (int i=0; i<traits::nrow(x); i++) {
        f(xptr[i]);
      }
    }
    return x;
  }

  ////

  // template <typename T>
  // T& dplusa(T& x, double a) {
  //   return mapapply(x, [=](double& value){ value += a; });
  // }

  // template <typename VectorT>
  // VectorT& dmulti(VectorT& x, const VectorT& y) {
  //   return map_vector(x, [=](int i, double& value){ value *= y[i]; });
  // }
  //
  // template <typename VectorT>
  // VectorT& ddiv(VectorT& x, const VectorT& y) {
  //   return map_vector(x, [=](int i, double& value){ value /= y[i]; });
  // }
  //
  // template <typename MatrixT, typename VectorT>
  // MatrixT& dmultiMN(MatrixT& x, const VectorT& y) {
  //   return map_matrix(x, [=](int i, int, double& value){ value *= y[i]; });
  // }
  //
  // template <typename MatrixT, typename VectorT>
  // MatrixT& dmultiMT(MatrixT& x, const VectorT& y) {
  //   return map_matrix(x, [=](int, int j, double& value){ value *= y[j]; });
  // }
  //
  // template <typename MatrixT, typename VectorT>
  // MatrixT& ddivMN(MatrixT& x, const VectorT& y) {
  //   return map_matrix(x, [=](int i, int, double& value){ value /= y[i]; });
  // }
  //
  // template <typename MatrixT, typename VectorT>
  // MatrixT& ddivMT(MatrixT& x, const VectorT& y) {
  //   return map_matrix(x, [=](int, int j, double& value){ value /= y[j]; });
  // }
  //
  // //
  //
  // template <typename MatrixT>
  // inline
  // MatrixT& dunif(MatrixT& x, double ufact, double& qv) {
  //   double* xptr = &x[0];
  //   double maxv = 0;
  //   for (int i=0; i<nrow(x); i++, xptr += ld(x)+1) {
  //     double tmp = std::abs(*xptr);
  //     if (tmp > maxv) {
  //       maxv = tmp;
  //     }
  //   }
  //   qv = maxv * ufact;
  //   dscal(1/qv, x);
  //
  //   xptr = &x[0];
  //   for (int i=0; i<nrow(x); i++, xptr += ld(x)+1) {
  //     *xptr += 1;
  //   }
  // }
}
