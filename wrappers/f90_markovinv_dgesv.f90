!
! wappers for markovinv_dgesv_mod.f90
!

subroutine f90_markovinv_dgesv_dense(trans, n, nrhs, Q, ldq, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, ldq, ldx, ldy
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(in) :: x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  ! TODO: check a cause of this error
  ! to avoid an error caused by a bug in gfortran?
  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_dense(trans, n, nrhs, Q, ldq, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_dense

subroutine f90_markovinv_dgesv_csr(trans, n, nrhs, spQ, rowptr, colind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(in) :: x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_csr(trans, n, nrhs, spQ, rowptr, colind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_csr

subroutine f90_markovinv_dgesv_csc(trans, n, nrhs, spQ, colptr, rowind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_csc(trans, n, nrhs, spQ, colptr, rowind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_csc

subroutine f90_markovinv_dgesv_coo(trans, n, nrhs, spQ, rowind, colind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_coo(trans, n, nrhs, spQ, rowind, colind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_coo
