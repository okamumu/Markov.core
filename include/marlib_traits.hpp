
namespace marlib {

  struct int_vector_tag{};
  struct double_vector_tag{};
  struct double_dense_matrix_tag{};
  struct double_csr_matrix_tag{};
  struct constant_value_tag{};

  struct eye_matrix{};

  template <class T>
  struct vector_traits {
    using value_type = typename T::value_type;
    static const value_type* begin(const T& v) { return &v[0]; }
    static const value_type* end(const T& v) { return &v[0] + v.size(); }
    static value_type* begin(T& v) { return &v[0]; }
    static value_type* end(T& v) { return &v[0] + v.size(); }
    static int size(const T& v) { return v.size(); }
    static int inc(const T&) { return 1; }
  };

  template <class T>
  struct dense_matrix_traits {
    using value_type = typename T::value_type;
    static const value_type* begin(const T& m) { return &m[0]; }
    static const value_type* end(const T& m) { return &m[0] + m.size(); }
    static value_type* begin(T& m) { return &m[0]; }
    static value_type* end(T& m) { return &m[0] + m.size(); }
    static int size(const T& m) { return m.size(); }
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int ld(const T& m) { return m.nrow(); }
  };

  template <class T>
  struct csr_matrix_traits {
    using value_type = typename T::value_type;
    static const value_type* begin(const T& m) { return &m[0]; }
    static const value_type* end(const T& m) { return &m[0] + m.size(); }
    static value_type* begin(T& m) { return &m[0]; }
    static value_type* end(T& m) { return &m[0] + m.size(); }
    static int size(const T& m) { return m.size(); }
    static int nnz(const T& m) { return m.nnz(); }
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static const int* rowptr(const T& m) { return m.rowptr(); }
    static const int* colind(const T& m) { return m.colind(); }
    static int* rowptr(T& m) { return m.rowptr(); }
    static int* colind(T& m) { return m.colind(); }
    static int origin(T& m) { return m.origin(); }
  };

  template <class T, class L>
  struct get_category;

  struct blas_level1_tag{};
  struct blas_level2_tag{};
  struct blas_level3_tag{};

  template <class L>
  struct get_category<int,L> {
    using type = constant_value_tag;
  };

  template <class L>
  struct get_category<double,L> {
    using type = constant_value_tag;
  };

  template <typename Value1, typename Value2>
  Value2& dcopy(const Value1& x, Value2& y) {
    return dcopy_impl(x, y,
      typename get_category<Value1,blas_level1_tag>::type(),
      typename get_category<Value2,blas_level1_tag>::type());
  }

  // const_value to ...

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_vector_tag) {
    using traits2 = vector_traits<Value2>;
    dblas::dfill(traits2::size(y), traits2::begin(y), traits2::inc(y), static_cast<typename traits2::value_type>(x));
    return y;
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_dense_matrix_tag) {
    using traits2 = dense_matrix_traits<Value2>;
    typename traits2::value_type* yptr = traits2::begin(y);
    for (int j=0; j<traits2::ncol(y); j++, yptr+=traits2::ld(y)) {
      dblas::dfill(traits2::nrow(y), yptr, 1, static_cast<typename traits2::value_type>(x));
    }
    return y;
  }

  // vector to ...

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_vector_tag) {
    using traits1 = vector_traits<Value1>;
    using traits2 = vector_traits<Value2>;
    assert(traits1::size(x) == traits2::size(y));
    dblas::dcopy(traits1::size(x), traits1::begin(x), traits1::inc(x), traits2::begin(y), traits2::inc(y));
    return y;
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_dense_matrix_tag) {
    using traits1 = vector_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    typename traits1::value_type* xptr = traits1::begin(x);
    typename traits2::value_type* yptr = traits2::begin(y);
    for (int j=0; j<traits2::ncol(y); j++, xptr+=traits1::inc(x)*traits2::nrow(y), yptr+=traits2::ld(y)) {
      dblas::dcopy(traits2::nrow(y), xptr, traits1::inc(x), yptr, 1);
    }
    return y;
  }

  // template <typename Value1, typename Value2>
  // Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_dense_matrix_tag) {
  //   using traits1 = dense_matrix_traits<Value1>;
  //   using traits2 = dense_matrix_traits<Value2>;
  //   dblas::dcopy(traits1::size(x), traits1::begin(x), traits1::inc(x), traits2::begin(y), traits2::inc(y));
  //   return y;
  // }
}
