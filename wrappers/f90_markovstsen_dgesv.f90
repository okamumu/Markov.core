!
! wrappers for markovstsen_dgesv_mod.f90
!

subroutine f90_markovstsen_dgesv_dense(n, Q, ldq, dQ, lddq, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, ldq, lddq, incp, incs
  double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_dense(n, Q, ldq, dQ, lddq, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_dense

subroutine f90_markovstsen_dgesv_csr(n, spQ, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
  integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_csr(n, spQ, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_csr

subroutine f90_markovstsen_dgesv_csc(n, spQ, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
  integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_csc(n, spQ, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_csc

subroutine f90_markovstsen_dgesv_coo(n, spQ, rowind0, colind0, nnz0, &
  dQ, rowind1, colind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
  integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_coo(n, spQ, rowind0, colind0, nnz0, &
  dQ, rowind1, colind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_coo
