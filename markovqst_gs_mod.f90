
module markovqst_gs
  implicit none
  private
  public markovqst_gs_dense
  public markovqst_gs_csr
  public markovqst_gs_csc
!   public markovqst_gs_coo

contains

  ! Description: quasi-stationary vector
  !
  !        x * Q = -gam * x,   x*1 = 1
  !
  !
  ! Q: lossy CTMC kernel
  ! xi = -Q*1: exit vector

  subroutine markovqst_gs_dense(n, Q, ldq, xi, incxi, &
    xstart, incxs, x, incx, &
    gam, maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    integer, intent(in) :: n, ldq, incxi, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
    double precision, intent(in) :: Q(1:ldq,1:n), xi(1:incxi,1:n), xstart(1:incxs,1:n)
    double precision, intent(out) :: x(1:incx,1:n), gam, rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: ddot
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    prevx(1:n) = 0.0d0
    bzero(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        gam = ddot(n, x, incx, xi, incxi)
        call gsstep_fwd_shift_dense('T', n, 1.0d0, Q, ldq, -gam, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0d0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx(1:n)/x(1,1:n)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovqst_gs_dense

  subroutine markovqst_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
    xstart, incxs, x, incx, &
    gam, maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    integer, intent(in) :: n, nnz, incxi, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
    double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n), xstart(1:incxs,1:n)
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n), gam, rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: ddot
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    prevx(1:n) = 0.0d0
    bzero(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        gam = ddot(n, x, incx, xi, incxi)
        call gsstep_fwd_shift_csr('T', n, 1.0d0, Q, rowptr, colind, nnz, -gam, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0d0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx(1:n)/x(1,1:n)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovqst_gs_csr

  subroutine markovqst_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
    xstart, incxs, x, incx, &
    gam, maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    integer, intent(in) :: n, nnz, incxi, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
    double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n), xstart(1:incxs,1:n)
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n), gam, rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: ddot
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    prevx(1:n) = 0.0d0
    bzero(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        gam = ddot(n, x, incx, xi, incxi)
        call gsstep_fwd_shift_csc('T', n, 1.0d0, Q, colptr, rowind, nnz, -gam, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0d0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx(1:n)/x(1,1:n)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovqst_gs_csc

end module markovqst_gs
