!
! wrappers for markovst_dgesv_mod.f90
!

subroutine f90_markovst_dgesv_dense(n, Q, ldq, x, incx, info)
  use markovst_dgesv
  integer, intent(in) :: n, ldq, incx
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(out) :: x(1:incx,1:n)
  integer, intent(out) :: info

  call markovst_dgesv_dense(n, Q, ldq, x, incx, info)
end subroutine f90_markovst_dgesv_dense

subroutine f90_markovst_dgesv_csr(n, spQ, rowptr, colind, nnz, x, incx, info)
  use markovst_dgesv
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)
  integer, intent(out) :: info

  call markovst_dgesv_csr(n, spQ, rowptr, colind, nnz, x, incx, info)
end subroutine f90_markovst_dgesv_csr

subroutine f90_markovst_dgesv_csc(n, spQ, colptr, rowind, nnz, x, incx, info)
  use markovst_dgesv
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)
  integer, intent(out) :: info

  call markovst_dgesv_csc(n, spQ, colptr, rowind, nnz, x, incx, info)
end subroutine f90_markovst_dgesv_csc

subroutine f90_markovst_dgesv_coo(n, spQ, rowind, colind, nnz, x, incx, info)
  use markovst_dgesv
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)
  integer, intent(out) :: info

  call markovst_dgesv_coo(n, spQ, rowind, colind, nnz, x, incx, info)
end subroutine f90_markovst_dgesv_coo
