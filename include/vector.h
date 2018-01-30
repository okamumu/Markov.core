/*
  dense_matrix
*/

namespace marlib {

  template <typename VectorT>
  class vector {

  using iterator = double*;

  public:
    vector(int size, VectorT v)
    : m_size(size), m_value(v) {}
    vector(const vector<VectorT>& v)
    : m_size(v.m_size), m_value(v.m_value) {}

    ~vector() {}

  private:
    int m_size;
    VectorT m_value;

  public:
    double& operator[](int i) { return m_value[i]; };
    const double& operator[](int i) const { return m_value[i]; };

    iterator begin() { return &m_value[0]; }
    iterator end() { return &m_value[0] + m_size; }
    const iterator begin() const { return m_value; }
    const iterator end() const { return m_value + m_size; }

    int size() const { return m_size; }

    ////// print
    std::ostream& print(std::ostream& os) const {
      for (int i=0; i<m_size; i++) {
        os << m_value[i] << " ";
      }
      os << std::endl;
      return os;
    }

    template <typename VectorTT>
    friend std::ostream& operator<< (std::ostream& os, const vector<VectorTT>& m);
  };

  template <typename VectorT>
  std::ostream& operator<<(std::ostream& os, const vector<VectorT>& m) {
    return m.print(os);
  }

}
