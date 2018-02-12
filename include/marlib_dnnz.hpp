namespace marlib {

  template <typename Value1>
  int dnnz(const Value1& x) {
    return dnnz_impl(x, typename get_category<Value1>::type());
  }

  // dense_matrix

  template <typename Value1>
  int dnnz_impl(const Value1& x, double_dense_matrix_tag) {
    using traits1 = dense_matrix_traits<Value1>;
    return dblas::dnnz(traits1::nrow(x), traits1::ncol(x), traits1::begin(x), traits1::ld(x));
  }

  // csr_matrix

  template <typename Value1>
  int dnnz_impl(const Value1& x, double_csr_matrix_tag) {
    using traits1 = csr_matrix_traits<Value1>;
    return traits1::nnz(x);
  }

}
