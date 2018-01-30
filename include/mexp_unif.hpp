/*
! Description: matrix exp with uniformization
!
!        y = exp(Q*t)
!
!        Q is uniformized to P and qv
!        t is involved in the Poisson probability vector.
*/

namespace marlib {

  template <typename MatrixT, typename VectorT, typename Poisson>
  VectorT& mexp_unif_trans(const MatrixT& P, double qv, const Poisson& pois, const VectorT& x, VectorT& y, double atol) {
    VectorT xi = clone(x);
    VectorT tmp = clone(x);

    dfill(y, 0);
    daxpy(pois(pois.left()), xi, y);
    for (int k=pois.left()+1; k<=pois.right(); k++) {
      dcopy(xi, tmp);
      dgemmTN(1, P, tmp, 0, xi);
      daxpy(pois(k), xi, y);
      if (dasum(xi) < atol)
        break;
    }
    dscal(1.0/pois.weight(), y);
    return y;
  }

  template <typename MatrixT, typename VectorT, typename Poisson>
  VectorT& mexp_unif_notrans(const MatrixT& P, double qv, const Poisson& pois, const VectorT& x, VectorT& y, double atol) {
    VectorT xi = clone(x);
    VectorT tmp = clone(x);

    dfill(y, 0);
    daxpy(pois(pois.left()), xi, y);
    for (int k=pois.left()+1; k<=pois.right(); k++) {
      dcopy(xi, tmp);
      dgemmNN(1, P, tmp, 0, xi);
      daxpy(pois(k), xi, y);
      if (dasum(xi) < atol)
        break;
    }
    dscal(1.0/pois.weight(), y);
    return y;
  }

}
