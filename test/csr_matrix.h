/*
  vector class
 */

namespace marlib {

  class csr_matrix {

  public:
    csr_matrix(int row_, int col_, int* rowptr_, int* colind_, double* value_, int nnz_, int origin_)
    : m_row(row_), m_col(col_), m_nnz(nnz_), m_origin(origin_), rowptr(rowptr_), colind(colind_), value(value_) {}

    ~csr_matrix() {}

  public:
    int m_row;
    int m_col;
    int m_nnz;
    int m_origin;
    int* rowptr;
    int* colind;
    double* value;

  public:
    double& operator[](int i) { return value[i]; };
    const double& operator[](int i) const { return value[i]; };

    int size() const { return m_nnz; }
    int nnz() const { return m_nnz; }
    int nrow() const { return m_row; }
    int ncol() const { return m_col; }
    int origin() const { return m_origin; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_row; i++) {
        for (int z=rowptr[i]-m_origin; z<rowptr[i+1]-m_origin; z++) {
          os << "(" << i + m_origin << "," << colind[z] + m_origin << ")=" << value[z] << std::endl;
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
}
