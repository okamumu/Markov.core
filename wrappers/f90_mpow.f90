!
! Wrappers for mpow_mod.f90
!

subroutine f90_mpow_dense(trans, n, MA, lda, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, lde, m
  double precision, intent(in) :: MA(1:lda,1:n)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_dense(trans, n, MA, lda, ME, lde, m, info)
end subroutine f90_mpow_dense

subroutine f90_mpow_csr(trans, n, spMA, rowptr, colind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_csr(trans, n, spMA, rowptr, colind, nnz, ME, lde, m, info)
end subroutine f90_mpow_csr

subroutine f90_mpow_csc(trans, n, spMA, colptr, rowind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_csc(trans, n, spMA, colptr, rowind, nnz, ME, lde, m, info)
end subroutine f90_mpow_csc

subroutine f90_mpow_coo(trans, n, spMA, rowind, colind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_coo(trans, n, spMA, rowind, colind, nnz, ME, lde, m, info)
end subroutine f90_mpow_coo
