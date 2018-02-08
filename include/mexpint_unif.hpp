/*
  ! Description: integral operation for matrix exp form;
  !
  !                    |t
  !        cME = cME + | exp(Q*s) ds
  !                    |0
  !
  !        ME = exp(Q*t)
  !
  !        Q is uniformized to P and qv
  !        t is involved in the Poisson probability vector.
  !        return value is ME
 */

namespace marlib {
  template <typename MatrixT, typename VectorT, typename Poisson>
  VectorT& mexpint_unif(char trans, const MatrixT& P, double qv,
    const Poisson& pois, const VectorT& x, VectorT& y, VectorT& cy,
    VectorT& xi, VectorT& tmp) {

    int left = pois.left();
    int right = pois.right();
    std::unique_ptr<double[]> ptr_cpoi(new double [right - left + 1]);
    double* cpoi = ptr_cpoi.get() - left; // offset -> left
    cpoi[right] = 0;
    for (int k=right-1; k>=left; k--) {
      cpoi[k] = cpoi[k+1] + pois[k];
    }

    dcopy(x, xi);
    dfill(y, 0.0);
    dscal(qv * pois.weight(), cy);
    daxpy(pois[left], xi, y);
    daxpy(cpoi[left], xi, cy);
    for (int k=left+1; k<=right; k++) {
      dcopy(xi, tmp);
      dgemm(trans, 'N', 1.0, P, tmp, 0.0, xi);
      daxpy(pois[k], xi, y);
      daxpy(cpoi[k], xi, cy);
    }
    dscal(1/pois.weight(), y);
    dscal(1/(qv * pois.weight()), cy);
    return y;
  }

}
