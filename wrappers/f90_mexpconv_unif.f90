!
! wrappers for mexpconv_unif_mod.f90
!

subroutine f90_mexpconv_unif_dense_vec(transQ, transH, n, P, ldp, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, H, ldh, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, ldp, left, right, incx, incy, incz, ldh
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), H(1:ldh,1:n)

  call mexpconv_unif_dense_vec(transQ, transH, n, P, ldp, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, H, ldh, atol)
end subroutine f90_mexpconv_unif_dense_vec

subroutine f90_mexpconv_unif_csr_vec(transQ, transH, n, spP, rowptr, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_csr_vec(transQ, transH, n, spP, rowptr, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_csr_vec

subroutine f90_mexpconv_unif_csc_vec(transQ, transH, n, spP, colptr, rowind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_csc_vec(transQ, transH, n, spP, colptr, rowind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_csc_vec

subroutine f90_mexpconv_unif_coo_vec(transQ, transH, n, spP, rowind, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_coo_vec(transQ, transH, n, spP, rowind, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_coo_vec
