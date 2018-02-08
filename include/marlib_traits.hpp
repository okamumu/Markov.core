
namespace marlib {

  struct int_vector_type{};
  struct double_vector_type{};
  struct double_dense_matrix_type{};
  struct double_dense_matrix_simple_type{}; // ld == nrow
  struct double_csr_matrix_type{};

  template <class T, typename ValueT>
  struct vector_traits {
    using value_type = ValueT;
    static const value_type* begin(const T& v) { return &v[0]; }
    static const value_type* end(const T& v) { return &v[0] + v.size(); }
    static value_type* begin(T& v) { return &v[0]; }
    static value_type* end(T& v) { return &v[0] + v.size(); }
    static int size(const T& v) { return v.size(); }
    static int inc(const T& v) { return 1; }
  }

  template <class T, typename ValueT>
  struct dense_matrix_traits {
    using value_type = ValueT;
    static const value_type* begin(const T& v) { return &v[0]; }
    static const value_type* end(const T& v) { return &v[0] + v.size(); }
    static value_type* begin(T& m) { return &m[0]; }
    static value_type* end(T& m) { return &m[0] + v.size(); }
    static int size(const T& m) { return m.size(); }
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int ld(const T& m) { return m.nrow(); }
  }

  template <class T, typename ValueT>
  struct csr_matrix_traits {
    using value_type = ValueT;
    static const value_type* begin(const T& v) { return &v[0]; }
    static const value_type* end(const T& v) { return &v[0] + v.size(); }
    static value_type* begin(T& m) { return &m[0]; }
    static value_type* end(T& m) { return &m[0] + v.size(); }
    static int size(const T& m) { return m.size(); }
    static int nnz(const T& m) { return m.nnz(); }
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static const int* rowptr(const T& m) { return m.rowptr(); }
    static const int* colind(const T& m) { return m.colind(); }
    static int* rowptr(T& m) { return m.rowptr(); }
    static int* colind(T& m) { return m.colind(); }
    static int origin(T& m) { return m.origin(); }
  }

  template <class T>
  struct get_category;

  template <typename Value1, typename Value2>
  Value2& dcopy(const Vector1& x, Value2& y) {
    return dcopy_impl(x, y,
      typename get_category<Value1>::type(),
      typename get_category<Value2>::type());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_type, double_vector_type) {
    using traits1 = vector_traits<Value1>;
    using traits2 = vector_traits<Value2>;
    dblas::dcopy(traits1::size(x), traits1::begin(x), traits1::inc(x), traits2::begin(y), traits2::inc(y));
    return y;
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_type, double_dense_matrix_type) {
    using traits1 = dense_matrix_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    dblas::dcopy(traits1::size(x), traits1::begin(x), traits1::inc(x), traits2::begin(y), traits2::inc(y));
    return y;
  }
}
