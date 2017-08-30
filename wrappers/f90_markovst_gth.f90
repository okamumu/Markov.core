!
! wrappers for markovst_gth_mod.f90
!

subroutine f90_markovst_gth_dense(n, Q, ldq, x, incx)
  use markovst_gth
  integer, intent(in) :: n, ldq, incx
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_dense(n, Q, ldq, x, incx)
end subroutine f90_markovst_gth_dense

subroutine f90_markovst_gth_csr(n, spQ, rowptr, colind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_csr(n, spQ, rowptr, colind, nnz, x, incx)
end subroutine f90_markovst_gth_csr

subroutine f90_markovst_gth_csc(n, spQ, colptr, rowind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_csc(n, spQ, colptr, rowind, nnz, x, incx)
end subroutine f90_markovst_gth_csc

subroutine f90_markovst_gth_coo(n, spQ, rowind, colind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_coo(n, spQ, rowind, colind, nnz, x, incx)
end subroutine f90_markovst_gth_coo
