! Poisson Probability Computation
!   Hiroyuki Okamura, Hiroshima University
!   okamu@rel.hiroshima-u.ac.jp

module poisson
  implicit none
!  double precision, private, parameter :: log2piOver2 = 0.9189385332046727417803297364056176398d0
  double precision, private, parameter :: log2piOver2 = log(2.0d0 * atan(1.0d0) * 4.0d0) / 2.0d0
  double precision, private, parameter :: biselection_tol = 1.0d-8
  double precision, private, parameter :: PUpperBound = 1.0d-4
  double precision, private, parameter :: PLowerBound = 1.0d-200
  double precision, private, parameter :: logPUpperBound = log(PUpperBound)
  double precision, private, parameter :: logPLowerBound = log(PLowerBound)
  double precision, private, parameter :: QLowerBound = 3.0d0 ! floor(-sqrt(2) * ierfc(2*(1-PUpperBound)))
  double precision, private, parameter :: QUpperBound = 31.0d0 ! ceiling(-sqrt(2) * ierfc(2*(1-PLowerBound)))
  double precision, private, parameter :: minlambda = 3.0d0 ! threshold to change rightbound
  integer, private, parameter :: rightmax = 23 ! sum_{k=0}^rightmax poi_pdf(minlambda, k) < 1.0 - 1.0d-13
  private normalt, normalq
  public poisson_rightbound, poisson_prob

contains

  ! Description
  !   return a tail probability of standard normal distribution
  ! Parameters
  !   q: input value (in)

  function normalt(q) result(result)

    double precision, intent(in) :: q
    double precision :: q2, tmp, sum
    double precision :: result

    q2 = q*q
    tmp = q
    sum = 1.0d0 / tmp
    tmp = tmp * q2
    sum = sum - 1.0d0 / tmp
    tmp = tmp * q2
    sum = sum + 3.0d0 / tmp
    tmp = tmp * q2
    sum = sum - 15.0d0 / tmp
    tmp = tmp * q2
    sum = sum + 105.0d0 / tmp
    result = log(sum) - q2/2.0d0 - log2piOver2
  end function normalt

  ! Description
  !   return a quantile of standard normal distribution
  ! Parameters
  !   p: probability for the quantile (in)
  !   tol: tolerance error (in)
  !   info: computation status (out)
  !           0: success
  !          -1: p is not appropriate
  !

  function normalq(p) result(q)
    double precision, intent(in) :: p
    double precision :: q
    double precision :: lp, ql, qu, fm

    lp = log(p)
    if ( (lp > logPUpperBound) .or. (lp < logPLowerBound) ) then
      call exit(1)
    end if

    ql = QLowerBound
    qu = QUpperBound
    q = (ql + qu) / 2.0d0
    fm = normalt(q) - lp
    do while ( abs(fm) > biselection_tol )
      if (fm > 0) then
        ql = q
      else
        qu = q
      end if
      q = (ql + qu)/2.0d0
      fm = normalt(q) - lp
    end do
  end function normalq

  ! Description:
  !   return the right bound of Poisson range
  !              for a given tolerance error
  ! Parameters:
  !   lambda: Poisson rate (in)
  !   eps: tolerance error (in)
  !   info: computation status (out)

  function poisson_rightbound(lambda, eps) result(right)
    double precision, intent(in) :: lambda, eps
    integer :: right
    
    integer :: k
    double precision :: z, tmp
    double precision :: total

    if (lambda < minlambda) then
      tmp = exp(-lambda)
      total = tmp
      right = 0
      do k = 1, rightmax
        right = right + 1
        tmp = tmp * lambda / right
        total = total + tmp
        if (total + eps >= 1.0d0) exit
      end do
    else
      z = normalq(eps)
      right =  floor((z + sqrt(4.0d0 * lambda - 1.0d0))**2 / 4.0d0 + 1.0d0)
    end if
  end function poisson_rightbound

  ! Description & parameters :
  !   poisson_prob(lambda, left, right, prob, weight)
  !     lambda, real*8: Poisson parameter
  !     left, integer: left bound
  !     right, integer: right bound
  !     prob, real*8: Poisson probabilities. The range must be left:right
  !     weight, real*8: weight value

  subroutine poisson_prob(lambda, left, right, prob, weight)

    integer, intent(in) :: left, right
    double precision, intent(in) :: lambda
    double precision, intent(out) :: prob(left:right), weight
    
    ! work variables
    integer :: mode, j, t, s

    mode = floor(lambda)
    if (mode >= 1) then
     prob(mode) = exp(-lambda + real(mode, kind=8) * log(lambda) &
      - log2piOver2 &
      - (real(mode, kind=8) + 1.0d0/2.0d0) &
      * log(real(mode, kind=8)) + real(mode, kind=8))
    else
      prob(mode) = exp(-lambda)
    end if
    ! -- down --
    do j = mode, left+1, -1
      prob(j-1) = (real(j, kind=8)/lambda)*prob(j)
    end do
    ! -- up --
    do j = mode, right-1
      prob(j+1) = (lambda/real(j+1, kind=8))*prob(j)
    end do
    ! -- compute W --
    weight = 0.0
    s = left
    t = right
    do while (s < t)
      if (prob(s) <= prob(t)) then
        weight = weight + prob(s)
        s = s + 1
      else
        weight = weight + prob(t)
        t = t - 1
      end if
    end do
    weight = weight + prob(s)
  end subroutine poisson_prob

end module poisson

