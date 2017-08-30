!
! wrappers for arnoldi_mod.f90
!

subroutine f90_arnoldi_dense(trans, n, A, lda, x, incx, m, H, ldh, V, ldv, &
  beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:lda,1:n), x(1:incx,1:n), tol
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_dense(trans, n, A, lda, x, incx, m, H, ldh, V, ldv, &
    beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_dense

subroutine f90_arnoldi_csr(trans, n, A, rowptr, colind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_csr(trans, n, A, rowptr, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_csr

subroutine f90_arnoldi_csc(trans, n, A, colptr, rowind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_csc(trans, n, A, colptr, rowind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_csc

subroutine f90_arnoldi_coo(trans, n, A, rowind, colind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_coo(trans, n, A, rowind, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_coo
