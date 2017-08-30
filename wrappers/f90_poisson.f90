!
! wrappers for poisson_mod.f90
!

subroutine f90_poisson_rightbound(lambda, eps, right)
  use poisson
  double precision, intent(in) :: lambda, eps
  integer, intent(out) :: right

  right = poisson_rightbound(lambda, eps)
end subroutine f90_poisson_rightbound

subroutine f90_poisson_prob(lambda, left, right, prob, weight)
  use poisson
  integer, intent(in) :: left, right
  double precision, intent(in) :: lambda
  double precision, intent(out) :: prob(left:right), weight

  call poisson_prob(lambda, left, right, prob, weight)
end subroutine f90_poisson_prob
