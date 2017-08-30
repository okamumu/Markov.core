
module erlang_dist
  implicit none
  integer, private, parameter :: erl_max = 10
  
contains

  function logpdf_erlang(n, rate, x) result(result)
    use gamma
    integer, intent(in) :: n
    double precision, intent(in) :: rate, x
    double precision :: result
    result = n * log(rate) + (n-1) * log(x) - rate * x - logfact(n-1)
  end function logpdf_erlang

  function pdf_erlang(n, rate, x) result(result)
    integer, intent(in) :: n
    double precision, intent(in) :: rate, x
    double precision :: result
    result = exp(logpdf_erlang(n, rate, x))
  end function pdf_erlang

  function cdf_erlang(n, rate, x) result(result)
    use gamma
    use gamma_dist
    integer, intent(in) :: n
    double precision, intent(in) :: rate, x
    double precision :: result
    integer :: i
    double precision :: w

    if (n <= erl_max) then
      w = 1.0d0
      do i = 1, n-1
        w = w + (rate* x)**i / fact(i)
      end do
      result = 1.0d0 - exp(-rate*x) * w
    else
      result = cdf_gamma(n+1.0d0, rate, x, loggamma(n+1.0d0))
    end if
  end function cdf_erlang

  function ccdf_erlang(n, rate, x) result(result)
    use gamma
    use gamma_dist
    integer, intent(in) :: n
    double precision, intent(in) :: rate, x
    double precision :: result
    integer :: i
    double precision :: w

    if (n <= erl_max) then
      w = 1.0d0
      do i = 1, n-1
        w = w + (rate* x)**i / fact(i)
      end do
      result = exp(-rate*x) * w
    else
      result = ccdf_gamma(n+1.0d0, rate, x, loggamma(n+1.0d0))
    end if
  end function ccdf_erlang

end module erlang_dist


