!
! wrappers for markovst_sor_mod.f90
!

subroutine f90_markovst_sor_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_dense

subroutine f90_markovst_sor_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_csr

subroutine f90_markovst_sor_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_csc
