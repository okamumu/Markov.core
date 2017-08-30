!
! wrappers for unif_matrix_mod.f90
!

subroutine f90_unif_dense(n, Q, ldq, P, ldp, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, ldq, ldp
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(out) :: P(1:ldp,1:n)

  call unif_dense(n, Q, ldq, P, ldp, qv, ufact)
end subroutine f90_unif_dense

subroutine f90_unif_csr(n, spQ, rowptr, colind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: rowptr(1:*), colind(1:*)

  call unif_csr(n, spQ, rowptr, colind, nnz, spP, qv, ufact)
end subroutine f90_unif_csr

subroutine f90_unif_csc(n, spQ, colptr, rowind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: colptr(1:*), rowind(1:*)

  call unif_csc(n, spQ, colptr, rowind, nnz, spP, qv, ufact)
end subroutine f90_unif_csc

subroutine f90_unif_coo(n, spQ, rowind, colind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: rowind(1:*), colind(1:*)

  call unif_coo(n, spQ, rowind, colind, nnz, spP, qv, ufact)
end subroutine f90_unif_coo
