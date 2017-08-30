
module markovst_gs
  implicit none
  private
  public markovst_gs_dense
  public markovst_gs_csr
  public markovst_gs_csc
!   public markovst_gs_coo

contains

  subroutine markovst_gs_dense(n, Q, ldq, xstart, incxs, x, incx, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use spblas
    use gsstep
    integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
    double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
    double precision, intent(out) :: x(1:incx,1:n), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call gsstep_fwd_shift_dense('T', n, 1.0d0, Q, ldq, 0.0d0, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx / x(1,1:n)))

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
  end subroutine markovst_gs_dense

  subroutine markovst_gs_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
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

    integer :: k
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call gsstep_fwd_shift_csr('T', n, 1.0d0, Q, rowptr, colind, nnz, 0.0d0, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx / x(1,1:n)))

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
  end subroutine markovst_gs_csr

  subroutine markovst_gs_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
    double precision, intent(in) :: rtol
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

    integer :: k
    double precision :: tmp, prevx(1:n), bzero(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0

    iter = 0
    do
      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call gsstep_fwd_shift_csc('T', n, 1.0d0, Q, colptr, rowind, nnz, 0.0d0, 1.0d0, bzero, 1, x, incx)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)
      iter = iter + steps

      rerror = maxval(abs(prevx / x(1,1:n)))

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
  end subroutine markovst_gs_csc

end module markovst_gs
