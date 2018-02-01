/*
  dense_matrix
*/

namespace marlib {

  class dense_matrix {

  public:
    dense_matrix(int nrow, int ncol, double* v)
    : m_size(nrow*ncol), m_row(nrow), m_col(ncol), m_value(v) {}
    ~dense_matrix() {}

  private:
    int m_size;
    int m_row;
    int m_col;
    double* m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

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

    friend std::ostream& operator<< (std::ostream& os, const dense_matrix& m);
  };

  inline
  std::ostream& operator<<(std::ostream& os, const dense_matrix& m) {
    return m.print(os);
  }

}
