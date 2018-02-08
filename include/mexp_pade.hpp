
namespace marlib {

  template <typename MatrixT, typename DenseMatrixT>
  DenseMatrixT& mexp_pade(char trans, const MatrixT& MA, double t,
    DenseMatrixT& ME, double eps,
    DenseMatrixT& MD, DenseMatrixT& MN, DenseMatrixT& MX, DenseMatrixT& MT) {

    dcopy(MA, ME);
    dscal(t, ME);
    double norma = damax(ME);

    int j = static_cast<int>(log(norma) / log(2.0));
    j = (0 < 1+j) ? 1+j : 0;
    dscal(1/std::ldexp(1.0,j), ME);

    int q = 1;
    double tolerr = 1.0 / 6.0;
    while (tolerr > eps/norma) {
      tolerr /= 16.0 * (3.0 + 4.0 * q * (2.0 + q));
      q++;
    }

    double c = 1.0;
    dfill(MD, eye());
    dfill(MN, eye());
    dfill(MX, eye());
    dfill(MT, 0.0);

    for (int k=1; k<=q; k++) {
      c *= (q - k + 1.0) / ((2.0 * q - k + 1.0) * k);
      dgemm(trans, 'N', 1.0, ME, MX, 0.0, MT);
      dcopy(MT, MX);
      daxpy(c, MX, MN);
      if (k % 2 == 0) {
        daxpy(c, MX, MD);
      } else {
        daxpy(-c, MX, MD);
      }
    }
    ME = MN;
    dgesv(MD, ME);
    for (int k=1; k<=j; k++) {
      dgemm('N', 'N', 1.0, ME, ME, 0.0, MT);
      dcopy(MT, ME);
    }
    return ME;
  }

}
