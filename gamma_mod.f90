! Gamma Probability Computation
!   Hiroyuki Okamura, Hiroshima University
!   okamu@rel.hiroshima-u.ac.jp

module gamma
  implicit none
!  double precision, parameter :: npi = 3.14159265358979324d0
  double precision, private, parameter :: log2pi = 1.83787706640934548d0
  double precision, private, parameter :: logpi = 1.14472988584940017d0
  integer, private, parameter :: CN = 8
  double precision, private, parameter :: &
    B0 = 1.0d0, B1 = -1.0d0 / 2.0d0, &
    B2 = 1.0d0 / 6.0d0, B4 = -1.0d0 / 30.0d0, &
    B6 = 1.0d0 / 42.0d0, B8 = -1.0d0 / 30.0d0, &
    B10 = 5.0d0 / 66.0d0, B12 = -691.0d0 / 2730.0d0, &
    B14 = 7.0d0 / 6.0d0, B16 = -3617.0d0 / 510.0d0
  integer, private, parameter :: factmax = 20
  double precision, private, parameter :: nfact(0:factmax) = (/ &
    1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, &
    120.0d0, 720.0d0, 5040.0d0, 40320.0d0, 362880.0d0, 3628800.0d0, &
    39916800.0d0, 479001600.0d0, 6227020800.0d0, 87178291200.0d0, &
    1307674368000.0d0, 20922789888000.0d0, 355687428096000.0d0, &
    6402373705728000.0d0, 121645100408832000.0d0, &
    2432902008176640000.0d0 /)
  double precision, private, parameter :: lognfact(0:factmax) = (/ &
    0.0d0, 0.0d0, &
    0.6931471805599453d0, 1.791759469228055d0, &
    3.1780538303479458d0, 4.787491742782046d0, 6.579251212010101d0, &
    8.525161361065415d0, 10.60460290274525d0, 12.801827480081469d0, &
    15.104412573075516d0, 17.502307845873887d0, 19.987214495661885d0, &
    22.552163853123425d0, 25.19122118273868d0, 27.89927138384089d0, &
    30.671860106080672d0, 33.50507345013689d0, 36.39544520803305d0, &
    39.339884187199495d0, 42.335616460753485d0 /)
  public fact, logfact, loggamma, psi, polygamma

contains
  function fact(s)
    integer, intent(in) :: s
    double precision :: fact
    if (s <= factmax) then
       fact = nfact(s)
    else
       fact = exp(loggamma(dble(s + 1)))
    end if
  end function fact

  function logfact(s)
    integer, intent(in) :: s
    double precision :: logfact
    if (s <= factmax) then
       logfact = lognfact(s)
    else
       logfact = loggamma(dble(s + 1))
    end if
  end function logfact

  function loggamma(s)
    double precision, intent(in) :: s
    double precision :: loggamma
    double precision :: x, v, w
    x = s
    v = 1.0d0
    do while (x < CN)
       v = v * x
       x = x + 1.0d0
    end do
    w = 1.0d0 / (x * x)
    loggamma = ((((((((B16 / (16.0d0 * 15.0d0)) * w &
         + (B14 / (14.0d0 * 13.0d0))) * w &
         + (B12 / (12.0d0 * 11.0d0))) * w &
         + (B10 / (10.0d0 * 9.0d0))) * w &
         + (B8 / (8.0d0 * 7.0d0))) * w &
         + (B6 / (6.0d0 * 5.0d0))) * w &
         + (B4 / (4.0d0 * 3.0d0))) * w &
         + (B2 / (2.0d0 * 1.0d0))) / x &
         + 0.5d0 * log2pi - log(v) - x + (x - 0.5d0) * log(x)
    return
  end function loggamma

  function psi(s)
    double precision, intent(in) :: s
    double precision :: psi
    double precision :: x, v, w
    x = s
    v = 0.0d0
    do while (x < CN)
       v = v + 1.0d0 / x
       x = x + 1.0d0
    end do
    w = 1.0d0 / (x * x)
    v = v + ((((((((B16 / 16.0d0) * w &
         + (B14 / 14.0d0)) * w &
         + (B12 / 12.0d0)) * w &
         + (B10 / 10.0d0)) * w &
         + (B8 / 8.0d0)) * w &
         + (B6 / 6.0d0)) * w &
         + (B4 / 4.0d0)) * w &
         + (B2 / 2.0d0)) * w &
         + 0.5d0 / x
    psi = log(x) - v
  end function psi

  function polygamma(n, s)
    integer, intent(in) :: n
    double precision, intent(in) :: s
    double precision :: polygamma
    integer :: k
    double precision :: x, t, u, v, w
    x = s
    u = 1.0d0
    do k = 1-n, -1
       u = u * k
    end do
    v = 0.0d0
    do while (x < CN)
       v = v + 1.0d0 / x**(n+1)
       x = x + 1.0d0
    end do
    w = x * x
    t = (((((((B16 &
         * (n + 15.0d0) * (n + 14.0d0) / (16.0d0 * 15.0d0 * w) + B14) &
         * (n + 13.0d0) * (n + 12.0d0) / (14.0d0 * 13.0d0 * w) + B12) &
         * (n + 11.0d0) * (n + 10.0d0) / (12.0d0 * 11.0d0 * w) + B10) &
         * (n +  9.0d0) * (n +  8.0d0) / (10.0d0 *  9.0d0 * w) + B8) &
         * (n +  7.0d0) * (n +  6.0d0) / ( 8.0d0 *  7.0d0 * w) + B6) &
         * (n +  5.0d0) * (n +  4.0d0) / ( 6.0d0 *  5.0d0 * w) + B4) &
         * (n +  3.0d0) * (n +  2.0d0) / ( 4.0d0 *  3.0d0 * w) + B2) &
         * (n +  1.0d0) * n / (2.0d0 * 1.0d0 * w) + 0.5d0 * n / x + 1.0d0
    polygamma = u * (t / x**n + n*v)
  end function polygamma

end module gamma
