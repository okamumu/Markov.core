
module gamma_dist
  use gamma
  implicit none
  private p_gamma, q_gamma
contains

  function logpdf_gamma(shape, rate, x, lgamshape) result(result)
    double precision, intent(in) :: shape, rate, x, lgamshape
    double precision :: result
    result = shape * log(rate) + (shape-1.0d0) * log(x) - rate * x - lgamshape
  end function logpdf_gamma

  function pdf_gamma(shape, rate, x, lgamshape) result(result)
    double precision, intent(in) :: shape, rate, x, lgamshape
    double precision :: result
    result = exp(logpdf_gamma(shape, rate, x, lgamshape))
  end function pdf_gamma

  function p_gamma(a, x, loggamma_a)
    double precision, intent(in) :: a, x, loggamma_a
    double precision :: p_gamma
    integer :: k
    double precision :: res, term, previous
    if (x >= 1.0d0 + a) then
       p_gamma = 1.0d0 - q_gamma(a, x, loggamma_a)
       return
    end if
    if (x == 0.0d0) then
       p_gamma = 0.0d0
       return
    end if
    term = exp(a * log(x) - x - loggamma_a) / a
    res = term
    do k = 1, 1000
       term = term * x / (a+k)
       previous = res
       res = res + term
       if (res == previous) then
          p_gamma = res
          return
       end if
    end do
    p_gamma = res
  end function p_gamma

  function q_gamma(a, x, loggamma_a)
    double precision, intent(in) :: a, x, loggamma_a
    double precision :: q_gamma
    integer :: k
    double precision :: res, w, temp, previous, la, lb
    la = 1.0d0
    lb = 1.0d0 + x - a
    if (x < 1.0d0 + a) then
       q_gamma = 1.0d0 - p_gamma(a, x, loggamma_a)
       return
    end if
    w = exp(a * log(x) - x - loggamma_a)
    res = w / lb
    do k = 2, 1000
       temp = ((k-1.0d0-a)*(lb-la) + (x+k)*lb)/k
       la = lb
       lb = temp
       w = w * (k-1.0d0-a)/k
       temp = w / (la * lb)
       previous = res
       res = res + temp
       if (res == previous) then
          q_gamma = res
          return
       end if
    end do
    q_gamma = res
  end function q_gamma

  function cdf_gamma(shape, rate, x, lgamshape) result(result)
    double precision, intent(in) :: shape, rate, x, lgamshape
    double precision :: result
    result = p_gamma(shape, rate*x, lgamshape)
  end function cdf_gamma

  function ccdf_gamma(shape, rate, x, lgamshape) result(result)
    double precision, intent(in) :: shape, rate, x, lgamshape
    double precision :: result
    result = q_gamma(shape, rate*x, lgamshape)
  end function ccdf_gamma

end module gamma_dist


