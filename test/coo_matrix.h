/*
  vector class
 */

namespace marlib {

  class coo_matrix {

  public:
    coo_matrix(int row_, int col_, int* rowind_, int* colind_, double* value_, int nnz_, int origin_)
    : m_row(row_), m_col(col_), m_nnz(nnz_), m_origin(origin_), rowind(rowind_), colind(colind_), value(value_) {}

    ~coo_matrix() {}

  public:
    int m_row;
    int m_col;
    int m_nnz;
    int m_origin;
    int* rowind;
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
      for (int z=0; z<m_nnz; z++) {
        os << "(" << rowind[z] << "," << colind[z] << ")=" << value[z] << std::endl;
      }
      return os;
    }

    friend std::ostream& operator<< (std::ostream& os, const coo_matrix& m);

  };

  inline
  std::ostream& operator<<(std::ostream& os, const coo_matrix& m) {
    return m.print(os);
  }
}
