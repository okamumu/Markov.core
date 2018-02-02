/*
  vector class
 */

namespace marlib {

  class csc_matrix {

  public:
    csc_matrix(int row_, int col_, int* colptr_, int* rowind_, double* value_, int nnz_, int origin_)
    : m_row(row_), m_col(col_), m_nnz(nnz_), m_origin(origin_), colptr(colptr_), rowind(rowind_), value(value_) {}

    ~csc_matrix() {}

  public:
    int m_row;
    int m_col;
    int m_nnz;
    int m_origin;
    int* colptr;
    int* rowind;
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
      for (int j=0; j<m_col; j++) {
        for (int z=colptr[j]-m_origin; z<colptr[j+1]-m_origin; z++) {
          os << "(" << rowind[z] << "," << j + m_origin << ")=" << value[z] << std::endl;
        }
      }
      return os;
    }

    friend std::ostream& operator<< (std::ostream& os, const csc_matrix& m);

  };

  inline
  std::ostream& operator<<(std::ostream& os, const csc_matrix& m) {
    return m.print(os);
  }
}
