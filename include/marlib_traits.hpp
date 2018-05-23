
namespace marlib {

  struct constant_value_tag{};
  struct int_vector_tag{};
  struct double_vector_tag{};
  struct double_dense_matrix_tag{};
  struct double_csr_matrix_tag{};

  struct eye_matrix{};

  template <class T>
  struct get_category;

  template <>
  struct get_category<int> {
    using type = constant_value_tag;
  };

  template <>
  struct get_category<double> {
    using type = constant_value_tag;
  };

  template <class T>
  struct to_vector_ {
    using type = std::false_type;
  };

  template<typename T>
  struct add_const {
    using type = T const;
  };

  template<typename T>
  struct add_reference {
    using type = T&;
  };

  template<typename T>
  struct add_const_reference {
    using type = const T&;
  };

  template<typename T>
  struct add_pointer {
    using type = T*;
  };

  template<typename T>
  struct add_const_pointer {
    using type = const T*;
  };

  template<bool Cond, typename Then, typename Else>
  struct if_ {
    using type = Then;
  };

  template<typename Then, typename Else>
  struct if_<false, Then, Else> {
    using type = Else;
  };

  template <class... Args>
  struct is_callable_impl {
    template <class F>
    static std::true_type
      check(decltype(std::declval<F>()(std::declval<Args>()...), (void)0)*);

    template <class F>
    static std::false_type check(...);
  };

  template <class F, class... Args>
  struct is_callable : decltype(is_callable_impl<Args...>::template check<F>(nullptr)) {};

  template <class T>
  struct base_traits {
    using value_type = typename T::value_type;
    using pointer_type = typename add_pointer<typename T::value_type>::type;
    using reference_type = typename add_reference<typename T::value_type>::type;
    using const_reference_type = typename add_const_reference<typename T::value_type>::type;
    using const_pointer_type = typename add_const_pointer<typename T::value_type>::type;
    static const_pointer_type begin(const T& v) { return &v[0]; }
    static const_pointer_type end(const T& v) { return &v[0] + v.size(); }
    static pointer_type begin(T& v) { return &v[0]; }
    static pointer_type end(T& v) { return &v[0] + v.size(); }
    static int size(const T& v) { return v.size(); }
  };

  template <class T>
  struct vector_traits : public base_traits<T> {
    // static int inc(const T& v) { return v.inc(); }
    static int inc(const T&) { return 1; }
  };

  template <class T>
  struct dense_matrix_traits : public base_traits<T> {
    using to_vector = typename to_vector_<T>::type;
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static int ld(const T& m) { return m.ld(); }
  };

  template <class T>
  struct csr_matrix_traits : public base_traits<T> {
    static int nnz(const T& m) { return m.nnz(); }
    static int nrow(const T& m) { return m.nrow(); }
    static int ncol(const T& m) { return m.ncol(); }
    static const int* rowptr(const T& m) { return m.rowptr(); }
    static const int* colind(const T& m) { return m.colind(); }
    static int* rowptr(T& m) { return m.rowptr(); }
    static int* colind(T& m) { return m.colind(); }
    static int origin(const T& m) { return m.origin(); }
  };

}
