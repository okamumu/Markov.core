namespace marlib {

  template <typename Value1, typename Value2>
  Value2& dcopy(const Value1& x, Value2& y) {
    return dcopy_impl(x, y,
      typename get_category<Value1>::type(),
      typename get_category<Value2>::type());
  }

  template <typename Value2>
  Value2& dcopy(eye_matrix, Value2& y) {
    return dcopy_eye_impl(y, typename get_category<Value2>::type());
  }

  // const_value to vector

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_vector_tag) {
    using traits2 = vector_traits<Value2>;
    dblas::dfill(traits2::size(y), traits2::begin(y), traits2::inc(y), static_cast<typename traits2::value_type>(x));
    return y;
  }

  // const_value to dense_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag t1, double_dense_matrix_tag t2) {
    using traits2 = dense_matrix_traits<Value2>;
    return dcopy_impl(x, y, t1, t2, typename traits2::to_vector());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_dense_matrix_tag, std::true_type) {
    std::cout << "call simple" << std::endl;
    return dcopy_impl(x, y, constant_value_tag(), double_vector_tag());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_dense_matrix_tag, std::false_type) {
    std::cout << "call not simple" << std::endl;
    using traits2 = dense_matrix_traits<Value2>;
    auto yptr = traits2::begin(y);
    for (int j=0; j<traits2::ncol(y); j++, yptr+=traits2::ld(y)) {
      dblas::dfill(traits2::nrow(y), yptr, 1, static_cast<typename traits2::value_type>(x));
    }
    return y;
  }

  // const_value to csr_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, constant_value_tag, double_csr_matrix_tag) {
    return dcopy_impl(x, y, constant_value_tag(), double_vector_tag());
  }

  // vector to vector

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_vector_tag) {
    using traits1 = vector_traits<Value1>;
    using traits2 = vector_traits<Value2>;
    assert(traits1::size(x) == traits2::size(y));
    dblas::dcopy(traits1::size(x), traits1::begin(x), traits1::inc(x), traits2::begin(y), traits2::inc(y));
    return y;
  }

  // vector to dense_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag t1, double_dense_matrix_tag t2) {
    using traits2 = dense_matrix_traits<Value2>;
    return dcopy_impl(x, y, t1, t2, typename traits2::to_vector());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_dense_matrix_tag, std::true_type) {
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_dense_matrix_tag, std::false_type) {
    using traits1 = vector_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    auto xptr = traits1::begin(x);
    auto yptr = traits2::begin(y);
    assert(traits1::size(x) == traits2::nrow(y)*traits2::ncol(y));
    for (int j=0; j<traits2::ncol(y); j++, xptr+=traits1::inc(x)*traits2::nrow(y), yptr+=traits2::ld(y)) {
      dblas::dcopy(traits2::nrow(y), xptr, traits1::inc(x), yptr, 1);
    }
    return y;
  }

  // vector to csr_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_vector_tag, double_csr_matrix_tag) {
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  // dense_matrix to vector

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag t1, double_vector_tag t2) {
    using traits1 = dense_matrix_traits<Value1>;
    return dcopy_impl(x, y, t1, t2, typename traits1::to_vector());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_vector_tag, std::true_type) {
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_vector_tag, std::false_type) {
    using traits1 = dense_matrix_traits<Value1>;
    using traits2 = vector_traits<Value2>;
    auto xptr = traits1::begin(x);
    auto yptr = traits2::begin(y);
    assert(traits2::size(y) == traits1::nrow(x)*traits1::ncol(x));
    for (int j=0; j<traits1::ncol(x); j++, xptr+=traits1::ld(x), yptr+=traits2::inc(y)*traits1::nrow(x)) {
      dblas::dcopy(traits1::nrow(x), xptr, 1, yptr, traits2::inc(y));
    }
    return y;
  }

  // dense_matrix to dense_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag t1, double_dense_matrix_tag t2) {
    std::cout << "call dcopy dmat dmat" << std::endl;
    using traits1 = dense_matrix_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    return dcopy_impl(x, y, t1, t2, typename traits1::to_vector(), typename traits2::to_vector());
  }

  template <typename Value1, typename Value2, typename BoolType1, typename BoolType2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_dense_matrix_tag, BoolType1, BoolType2) {
    std::cout << "call not simple mat" << std::endl;
    using traits1 = dense_matrix_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    auto xptr = traits1::begin(x);
    auto yptr = traits2::begin(y);
    assert(traits1::nrow(x) == traits2::nrow(y));
    assert(traits1::ncol(x) == traits2::ncol(y));
    for (int j=0; j<traits1::ncol(x); j++, xptr+=traits1::ld(x), yptr+=traits2::ld(y)) {
      dblas::dcopy(traits1::nrow(x), xptr, 1, yptr, 1);
    }
    return y;
  }

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_dense_matrix_tag, std::true_type, std::true_type) {
    std::cout << "call simple mat" << std::endl;
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  // dense_matrix to csr_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_dense_matrix_tag, double_csr_matrix_tag) {
    using traits1 = dense_matrix_traits<Value1>;
    using traits2 = csr_matrix_traits<Value2>;
    assert(traits1::nrow(x) == traits2::nrow(y));
    assert(traits1::ncol(x) == traits2::ncol(y));
    dblas::dense_to_csr(traits1::nrow(x), traits1::ncol(x), traits1::begin(x), traits1::ld(x),
      traits2::begin(y), traits2::rowptr(y), traits2::colind(y), traits2::nnz(y), traits2::origin(y));
    return y;
  }

  // csr_matrix to vector

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_csr_matrix_tag, double_vector_tag) {
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  // csr_matrix to dense_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_csr_matrix_tag, double_dense_matrix_tag) {
    using traits1 = csr_matrix_traits<Value1>;
    using traits2 = dense_matrix_traits<Value2>;
    assert(traits1::nrow(x) == traits2::nrow(y));
    assert(traits1::ncol(x) == traits2::ncol(y));
    dblas::csr_to_dense(traits1::nrow(x), traits1::ncol(x), traits1::begin(x), traits1::rowptr(x),
      traits1::colind(x), traits1::nnz(x), traits1::origin(x),
      traits2::begin(y), traits2::ld(y));
    return y;
  }

  // csr_matrix to csr_matrix

  template <typename Value1, typename Value2>
  Value2& dcopy_impl(const Value1& x, Value2& y, double_csr_matrix_tag, double_csr_matrix_tag) {
    return dcopy_impl(x, y, double_vector_tag(), double_vector_tag());
  }

  // eye for dense

  template <typename Value2>
  Value2& dcopy_eye_impl(Value2& y, double_dense_matrix_tag) {
    using traits2 = dense_matrix_traits<Value2>;
    auto yptr = traits2::begin(y);
    for (int j=0; j<traits2::ncol(y); j++, yptr+=traits2::ld(y)) {
      for (int i=0; i<traits2::nrow(y); i++) {
        if (i == j) {
          yptr[i] = 1;
        } else {
          yptr[i] = 0;
        }
      }
    }
    return y;
  }

  // eye for csr

  template <typename Value2>
  Value2& dcopy_eye_impl(Value2& y, double_csr_matrix_tag) {
    using traits2 = csr_matrix_traits<Value2>;
    auto value = traits2::begin(y);
    for (int i=0; i<traits2::nrow(y); i++) {
      for (int z=traits2::rowptr(y)[i]-traits2::origin(y); z<traits2::rowptr(y)[i+1]-traits2::origin(y); z++) {
        int j = traits2::colind(y)[z] - traits2::origin(y);
        std::cout << "(" << i << ", " << j << ")" << std::endl;
        if (i == j) {
          value[z] = 1;
        } else {
          value[z] = 0;
        }
      }
    }
    return y;
  }
}
