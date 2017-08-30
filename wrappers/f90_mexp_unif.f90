!
! wrappers for mexp_unif_mod.f90
!

subroutine f90_mexp_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, incx, incy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), x(1:*), atol
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_dense_vec

subroutine f90_mexp_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_csr_vec

subroutine f90_mexp_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_csc_vec

subroutine f90_mexp_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_coo_vec

subroutine f90_mexp_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, m, ldx, ldy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_dense_mat

subroutine f90_mexp_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_csr_mat

subroutine f90_mexp_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_csc_mat

subroutine f90_mexp_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_coo_mat
