/*
  vector class
 */

namespace marlib {

  template <typename VectorT, typename IntVectorT>
  class csr_matrix {

    using iterator = VectorT;
    using int_iterator = IntVectorT;

  public:
    csr_matrix(int row, int col,
      IntVectorT rowptr, IntVectorT colind, VectorT value, int nnz, int origin)
    : m_row(row), m_col(col), m_nnz(nnz), m_origin(origin), m_rowptr(rowptr), m_colind(colind), m_value(value) {}

    ~csr_matrix() {}

  private:
    int m_row;
    int m_col;
    int m_nnz;
    int m_origin;
    IntVectorT m_rowptr;
    IntVectorT m_colind;
    VectorT m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

    iterator begin() { return m_value; }
    iterator end() { return m_value + m_nnz; }
    const iterator begin() const { return m_value; }
    const iterator end() const { return m_value + m_nnz; }

    int size() const { return m_nnz; }
    int nnz() const { return m_nnz; }
    int nrow() const { return m_row; }
    int ncol() const { return m_col; }

    int origin() const { return m_origin; }
    int_iterator pbegin() { return m_rowptr; }
    int_iterator pend() { return m_rowptr + m_row + 1; };
    int_iterator cbegin() { return m_colind; }
    int_iterator cend() { return m_colind + m_nnz; }
    const int_iterator pbegin() const { return m_rowptr; }
    const int_iterator pend() const { return m_rowptr + m_row + 1; };
    const int_iterator cbegin() const { return m_colind; }
    const int_iterator cend() const { return m_colind + m_nnz; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_row; i++) {
        for (int z=m_rowptr[i]-m_origin; z<m_rowptr[i+1]-m_origin; z++) {
          os << "(" << i + m_origin << "," << m_colind[z] + m_origin << ")=" << m_value[z] << std::endl;
        }
      }
      return os;
    }

    template <typename VectorTT, typename IntVectorTT>
    friend std::ostream& operator<< (std::ostream& os, const csr_matrix<VectorTT,IntVectorTT>& m);

  };

  template <typename ValueT, typename RangeT>
  std::ostream& operator<<(std::ostream& os, const csr_matrix<ValueT,RangeT>& m) {
    return m.print(os);
  }

  template <typename VectorT, typename IntVectorT>
  inline
  csr_matrix<VectorT,IntVectorT>& dcopy(const dense_matrix<VectorT>& x, csr_matrix<VectorT,IntVectorT>& y) {
    auto r = y.pbegin();
    auto c = y.cbegin();
    dblas::dense_to_csr(x.nrow(), x.ncol(), &x[0], x.nrow(), &y[0], &r[0], &c[0], y.nnz(), y.origin());
    return y;
  }

  template <typename VectorT, typename IntVectorT>
  inline
  int nnz(const csr_matrix<VectorT,IntVectorT>& m) {
    return m.nnz();
  }
}
