!
! wrappers for gaussinte_mod.f90
!

subroutine f90_gaussinte_w(n, x, w, eps)
  use gaussinte
  integer, intent(in) :: n
  double precision, intent(out) :: x(1:n), w(1:n)
  double precision, intent(in) :: eps

  call gaussinte_w(n, x, w, eps)
end subroutine f90_gaussinte_w

function f90_gaussinte_fx(n, x, a, b, fx) result(c)
  use gaussinte
  integer, intent(in) :: n
  double precision, intent(in) :: x(1:n), a, b
  double precision, intent(out) :: fx(1:n)
  double precision :: c
  c = gaussinte_fx(n, x, a, b, fx)
end function f90_gaussinte_fx

function f90_gaussinte_fv(n, w, c, fv) result(s)
  use gaussinte
  integer, intent(in) :: n
  double precision, intent(in) :: w(1:n), c, fv(1:n)
  double precision :: s
  s = gaussinte_fv(n, w, c, fv)
end function f90_gaussinte_fv
