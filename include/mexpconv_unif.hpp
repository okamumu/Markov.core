/*
  ! Description: convolution integral operation for matrix exp form;
  !
  !                           |t
  ! transH(MH) = transH(MH) + | exp(transQ(Q)*s) * x * y' * exp(transQ(Q)*(t-s)) ds
  !                           |0
  !
  !        and
  !
  !        z = exp(transQ(Q)*t) * x
  !
  !        t is involved in the Poisson probability vector.
  !        qv is an uniformed parameter
  !        return value is z
 */

namespace marlib {

  template <typename MatrixT, typename VectorT, typename Poisson, typename VectorVectorT>
  VectorT& mexpconv_unif(char transQ, char transH,
    const MatrixT& P, double qv, const Poisson& pois,
    const VectorT& x, const VectorT& y, VectorT& z, MatrixT& H,
    VectorT& xi, VectorT& tmp, VectorVectorT& vc) {

    int left = pois.left();
    int right = pois.right();

    assert(left == 0 && right >= 1);

    const char ctranQ = (transQ == 'N' || transQ == 'n')?'T':'N';

    dfill(vc[right], 0.0);
    daxpy(pois(right), y, vc[right]);
    for (int l=right-1; l>=left+1; l--) {
      dgemm(ctranQ, 'N', 1.0, P, vc[l+1], 0.0, vc[l]);
      daxpy(pois[l], y, vc[l]);
    }

    dcopy(x, xi);
    dfill(z, 0.0);
    dscal(qv * pois.weight(), H);
    daxpy(pois[left], xi, z);
    dger(transH, 1.0, xi, vc[left+1], H);

    for (int l=left+1; l<=right-1; l++) {
      dcopy(xi, tmp);
      dgemm(transQ, 'N', 1.0, P, tmp, 0.0, xi);
      daxpy(pois[l], xi, z);
      dger(transH, 1.0, xi, vc[l+1], H);
    }
    dscal(1/pois.weight(), z);
    dscal(1/(qv * pois.weight()), H);
    return z;
  }
}
