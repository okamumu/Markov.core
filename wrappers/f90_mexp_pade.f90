!
! wrappers for mexp_pade_mod.f90
!

subroutine f90_mexp_pade_dense(trans, n, alpha, MA, lda, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, lde
  double precision, intent(in) :: alpha, MA(1:lda,1:n), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_dense(trans, n, alpha, MA, lda, ME, lde, eps)
end subroutine f90_mexp_pade_dense

subroutine f90_mexp_pade_csr(trans, n, alpha, spMA, rowptr, colind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_csr(trans, n, alpha, spMA, rowptr, colind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_csr

subroutine f90_mexp_pade_csc(trans, n, alpha, spMA, colptr, rowind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_csc(trans, n, alpha, spMA, colptr, rowind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_csc

subroutine f90_mexp_pade_coo(trans, n, alpha, spMA, rowind, colind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_coo(trans, n, alpha, spMA, rowind, colind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_coo
