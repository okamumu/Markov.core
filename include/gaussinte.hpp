
namespace marlib {

/*
!   Description: Gauss quadrature for the following integral

!      | b
!      |  f(x) dx
!      | a

!   gaussinte_w: make points and weights for n discrete points
!     n (in): the number of points. This is the size of both x and w.
!     x (out): x-axis points in the interval [-1, 1].
!     w (out): weights for the points.
!     eps (in): tolerance error.

!   gaussinte_fx: make x points for the interval [a,b]
!     n (in): the number of points.
!     x (in): the x points for interval [-1, 1].
!     a, b (in): lower and upper points for the integral
!     fx (out): x point for the interval [a, b]
!     return value: (b-a)/2

!   gaussinte_fv: compute the integral
!     n (in): the number of points.
!     w (in): weights for the x points in [-1,1]
!     c (in): (b-a)/2 ?
!     fv (in): function values at x points derived by gauss_inte_fx
!     return value: the interal value
*/

  template <typename VectorT>
  class gaussinte {
  public:
    gaussinte(int n, const VectorT& x, const VectorT& w, double eps);

    double comp_fx(double a, double b, VectorT fx) const;
    double comp_value(double c, VectorT fv) const;
    VectorT get_w() const;
    int get_n() const { return m_n; }

  private:
    int m_n;
    VectorT m_x;
    VectorT m_w;

    void comp_w(double eps);
  };

  /// methods

  template <typename VectorT>
  gaussinte<VectorT>::gaussinte(int n, const VectorT& x, const VectorT& w, double eps)
  : m_n(n), m_x(x), m_w(w) {
      comp_w(eps);
  }

  template <typename VectorT>
  void gaussinte<VectorT>::comp_w(double eps) {
    constexpr double pai = 3.14159265358979324;
    switch (m_n) {
      case 1:
        m_x[0] = 0.0;
        m_w[0] = 2.0;
        return;
      case 2:
        m_x[0] = -std::sqrt(1.0/3.0);
        m_w[0] = 1.0;
        m_x[1] = -m_x[1];
        m_w[1] = m_w[0];
        return;
      case 3:
        m_x[0] = -std::sqrt(0.6);
        m_w[0] = 5.0/9.0;
        m_x[1] = 0.0;
        m_w[1] = 8.0/9.0;
        m_x[2] = -m_x[0];
        m_w[2] = m_w[0];
        return;
    }
    int m = m_n / 2;
    double npai = pai / (m_n + 0.5);
    double p0, p1, p2, q0, q1, q2;
    for (int i=1; i<=m; i++) {
      double tmp = std::cos((i - 0.25) * npai);
      double dt = tmp;
      while (std::abs(dt) > std::abs(tmp) * eps) {
        p1 = tmp;
        p2 = (3.0 * tmp * tmp - 1.0) * 0.5;
        q1 = 1.0;
        q2 = 3.0 * tmp;
        for (int l=3; l<=m_n; l++) {
          p0 = p1;
          p1 = p2;
          p2 = ((l + l - 1) * tmp * p1 - (l-1) * p0) / l;
          q0 = q1;
          q1 = q2;
          q2 = ((l + l - 1) * (tmp * q1 + p1) - (l-1) * q0) / l;
        }
        dt = p2 / q2;
        tmp = tmp - dt;
      }
      m_x[i-1] = tmp;
      m_w[i-1] = 2.0 / (m_n * p1 * q2);
    }

    if (m_n % 2 == 1) {
      double tmp = m_n;
      for (int i=1; i<=m; i++) {
        tmp = tmp * (0.5 - i) / i;
      }
      m_x[m+1-1] = 0.0;
      m_w[m+1-1] = 2.0 / (tmp * tmp);
    }

    for (int i=1; i<=m; i++) {
      m_x[m_n+1-i-1] = m_x[i-1];
      m_x[i-1] = -m_x[i-1]; // reverse order
      m_w[m_n+1-i-1] = m_w[i-1];
    }
  }

  template <typename VectorT>
  VectorT gaussinte<VectorT>::get_w() const {
      return m_w;
  }

  template <typename VectorT>
  double gaussinte<VectorT>::comp_fx(double a, double b, VectorT fx) const {
    double t1 = (b - a)/2;
    double t2 = (b + a)/2;
    for (int i=0; i<m_n; i++) {
      fx[i] = t1 * m_x[i] + t2;
    }
    return t1;
  }

  template <typename VectorT>
  double gaussinte<VectorT>::comp_value(double c, VectorT fv) const {
    double s = 0;
    for (int i=0; i<m_n; i++) {
      s += m_w[i] * fv[i];
    }
    s *= c;
    return s;
  }
}
