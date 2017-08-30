!
! wrappers for markovstsen_gs_mod.f90
!

subroutine f90_markovstsen_gs_dense(n, Q, ldq, dQ, lddq, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, ldq, lddq, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n), pis(1:incp,1:n), sstart(1:incss,1:n)
  double precision, intent(out) :: s(1:incs,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovstsen_gs_dense(n, Q, ldq, dQ, lddq, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_dense

subroutine f90_markovstsen_gs_csr(n, Q, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz0), dQ(1:nnz1), pis(1:incp,1:n), sstart(1:incss,1:n)
  integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
  integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
  double precision, intent(out) :: s(1:incs,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovstsen_gs_csr(n, Q, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_csr

subroutine f90_markovstsen_gs_csc(n, Q, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz0), dQ(1:nnz1), pis(1:incp,1:n), sstart(1:incss,1:n)
  integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
  integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
  double precision, intent(out) :: s(1:incs,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovstsen_gs_csc(n, Q, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_csc
