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
  VectorT& mexp_unif(char trans, const MatrixT& P,
    double, const Poisson& pois, const VectorT& x, VectorT& y,
    VectorT& xi, VectorT& tmp, double atol) {

    dcopy(x, xi);
    dfill(y, 0.0);
    daxpy(pois(pois.left()), xi, y);
    for (int k=pois.left()+1; k<=pois.right(); k++) {
      dcopy(xi, tmp);
      dgemm(trans, 'N', 1.0, P, tmp, 0.0, xi);
      daxpy(pois[k], xi, y);
      if (dasum(xi) < atol)
        break;
    }
    dscal(1.0/pois.weight(), y);
    return y;
  }

}
