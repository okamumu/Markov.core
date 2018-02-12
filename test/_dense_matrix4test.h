/*
  dense_matrix
*/

namespace marlib {

  template <typename LD = std::true_type>
  class dense_matrix {
  public:
    using value_type = double;
    using ld_type = LD;

    dense_matrix(int nrow, int ncol, double* v, int ld = 0)
    : m_size(nrow*ncol), m_row(nrow), m_col(ncol), m_ld(ld), m_value(v) {
      if (m_ld == 0) {
        m_ld = nrow;
      }
    }

    ~dense_matrix() {}

  private:
    int m_size;
    int m_row;
    int m_col;
    int m_ld;
    double* m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

    int size() const { return m_size; }
    int nrow() const { return m_row; }
    int ncol() const { return m_col; }
    int ld() const { return m_ld; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_row; i++) {
        for (int j=0; j<m_col; j++) {
          os << m_value[i + m_ld * j] << " ";
        }
        os << std::endl;
      }
      return os;
    }

    template <typename T>
    friend std::ostream& operator<< (std::ostream& os, const dense_matrix<T>& m);
  };

  template <typename LD>
  inline
  std::ostream& operator<<(std::ostream& os, const dense_matrix<LD>& m) {
    return m.print(os);
  }

  // template <>
  // struct vector_traits<dense_matrix<std::true_type>> : public base_traits<dense_matrix<std::true_type>> {
  //   static int inc(const dense_matrix<std::true_type>&) { return 1; }
  // };
  //
  template <typename LD>
  struct get_category<dense_matrix<LD>> {
    using type = double_dense_matrix_tag;
  };

  template <>
  struct to_vector_<dense_matrix<std::true_type>> {
    using type = std::true_type;
  };

}
