/*
  vector class
 */

namespace marlib {

  class csr_matrix {

  public:
    using value_type = double;

    csr_matrix(int row_, int col_, int* rowptr_, int* colind_, double* value_, int nnz_, int origin_)
    : m_row(row_), m_col(col_), m_nnz(nnz_), m_origin(origin_), m_rowptr(rowptr_), m_colind(colind_), m_value(value_) {}

    ~csr_matrix() {}

  public:
    int m_row;
    int m_col;
    int m_nnz;
    int m_origin;
    int* m_rowptr;
    int* m_colind;
    double* m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

    int size() const { return m_nnz; }
    int nnz() const { return m_nnz; }
    int nrow() const { return m_row; }
    int ncol() const { return m_col; }
    int origin() const { return m_origin; }
    int* rowptr() { return m_rowptr; }
    const int* rowptr() const { return m_rowptr; }
    int* colind() { return m_colind; }
    const int* colind() const { return m_colind; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_row; i++) {
        for (int z=m_rowptr[i]-m_origin; z<m_rowptr[i+1]-m_origin; z++) {
          os << "(" << i + m_origin << "," << m_colind[z] << ")=" << m_value[z] << std::endl;
        }
      }
      return os;
    }

    friend std::ostream& operator<< (std::ostream& os, const csr_matrix& m);

  };

  inline
  std::ostream& operator<<(std::ostream& os, const csr_matrix& m) {
    return m.print(os);
  }

  template <>
  struct get_category<csr_matrix> {
    using type = double_csr_matrix_tag;
  };
}
