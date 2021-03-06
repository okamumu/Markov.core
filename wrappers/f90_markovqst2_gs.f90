!
! wrappers for markovqst2_gs_mod.f90
!

subroutine f90_markovqst2_gs_dense(n, Q, ldq, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, ldq, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_dense(n, Q, ldq, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_dense

subroutine f90_markovqst2_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, nnz, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_csr

subroutine f90_markovqst2_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, nnz, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_csc
