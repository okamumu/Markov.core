/*
  dense_matrix
*/

namespace marlib {

  template <typename VectorT>
  class dense_matrix {

  using iterator = double*;

  public:
    dense_matrix(int nrow, int ncol, VectorT v)
    : m_size(nrow*ncol), m_row(nrow), m_col(ncol), m_value(v) {}
    dense_matrix(const dense_matrix<VectorT>& m)
    : m_size(m.m_size), m_row(m.m_row), m_col(m.m_col), m_value(m.m_value) {}

    ~dense_matrix() {}

  private:
    int m_size;
    int m_row;
    int m_col;
    VectorT m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

    iterator begin() { return &m_value[0]; }
    iterator end() { return &m_value[0] + m_size; }
    const iterator begin() const { return m_value; }
    const iterator end() const { return m_value + m_size; }

    int size() const { return m_size; }
    int nrow() const { return m_row; }
    int ncol() const { return m_col; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_row; i++) {
        for (int j=0; j<m_col; j++) {
          os << m_value[i + m_row * j] << " ";
        }
        os << std::endl;
      }
      return os;
    }

    template <typename VectorTT>
    friend std::ostream& operator<< (std::ostream& os, const dense_matrix<VectorTT>& m);
  };

  template <typename VectorT>
  std::ostream& operator<<(std::ostream& os, const dense_matrix<VectorT>& m) {
    return m.print(os);
  }

}
