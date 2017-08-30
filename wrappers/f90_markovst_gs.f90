!
! wrappers for markovst_gs_mod.f90
!

subroutine f90_markovst_gs_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_dense

subroutine f90_markovst_gs_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_csr

subroutine f90_markovst_gs_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_csc
