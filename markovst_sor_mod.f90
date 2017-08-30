
module markovst_sor
  implicit none
  private
  public markovst_sor_dense
  public markovst_sor_csr
  public markovst_sor_csc
!   public markovst_gs_coo

contains

  subroutine markovst_sor_dense(n, Q, ldq, xstart, incxs, x, incx, &
    maxiter, rtol, steps0, iter, rerror, omega, info, callback, update_omega)
    use spblas
    use gsstep
    integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps0
    double precision, intent(in) :: rtol
    double precision, intent(inout) :: omega
    double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
    double precision, intent(out) :: x(1:incx,1:n), rerror
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

    integer :: k, steps
    double precision :: tmp, prevx(1:n), bzero(1:n)
    double precision :: prevx1(1:n), prevx2(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0
    prevx1(1:n) = 0.0d0
    prevx2(1:n) = 0.0d0

    iter = 0
    steps = steps0
    do
      ! callback
      call callback(iter, steps, rerror, omega, n, prevx, x, incx)

      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call dcopy(n, prevx1, 1, prevx2, 1)
        call dcopy(n, x, incx, prevx1, 1)
        call gsstep_fwd_shift_dense('T', n, 1.0d0, Q, ldq, 0.0d0, omega, bzero, 1, x, incx)
        call daxpy(n, -1.0d0, x, incx, prevx1, 1)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
        iter = iter + 1
        if (mod(iter, 10) == 0) then
          call dscal(n, tmp, prevx1, 1)
          call update_omega(iter, n, prevx1, prevx2, omega)
        end if
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)

      rerror = maxval(abs(prevx / x(1,1:n)))

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovst_sor_dense

  subroutine markovst_sor_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
    maxiter, rtol, steps0, iter, rerror, omega, info, callback, update_omega)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps0
    double precision, intent(in) :: rtol
    double precision, intent(inout) :: omega
    double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n), rerror
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

    integer :: k, steps
    double precision :: tmp, prevx(1:n), bzero(1:n)
    double precision :: prevx1(1:n), prevx2(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0
    prevx1(1:n) = 0.0d0
    prevx2(1:n) = 0.0d0

    iter = 0
    steps = steps0
    do
      ! callback
      call callback(iter, steps, rerror, omega, n, prevx, x, incx)

      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call dcopy(n, prevx1, 1, prevx2, 1)
        call dcopy(n, x, incx, prevx1, 1)
        call gsstep_fwd_shift_csr('T', n, 1.0d0, Q, rowptr, colind, nnz, 0.0d0, omega, bzero, 1, x, incx)
        call daxpy(n, -1.0d0, x, incx, prevx1, 1)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
        iter = iter + 1
        ! update
        if (mod(iter, 10) == 0) then
          call dscal(n, tmp, prevx1, 1)
          call update_omega(iter, n, prevx1, prevx2, omega)
        end if
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)

      rerror = maxval(abs(prevx / x(1,1:n)))

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovst_sor_csr

  subroutine markovst_sor_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
    maxiter, rtol, steps0, iter, rerror, omega, info, callback, update_omega)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps0
    double precision, intent(in) :: rtol
    double precision, intent(inout) :: omega
    double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n), rerror
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

    integer :: k, steps
    double precision :: tmp, prevx(1:n), bzero(1:n)
    double precision :: prevx1(1:n), prevx2(1:n)

    x(1,1:n) = xstart(1,1:n)
    bzero(1:n) = 0.0d0
    prevx(1:n) = 0.0d0
    prevx1(1:n) = 0.0d0
    prevx2(1:n) = 0.0d0

    iter = 0
    steps = steps0
    do
      ! callback
      call callback(iter, steps, rerror, omega, n, prevx, x, incx)

      call dcopy(n, x, incx, prevx, 1)
      do k = 1, steps
        call dcopy(n, prevx1, 1, prevx2, 1)
        call dcopy(n, x, incx, prevx1, 1)
        call gsstep_fwd_shift_csc('T', n, 1.0d0, Q, colptr, rowind, nnz, 0.0d0, omega, bzero, 1, x, incx)
        call daxpy(n, -1.0d0, x, incx, prevx1, 1)
        tmp = sum(x(1,1:n))
        call dscal(n, 1.0/tmp, x, incx)
        iter = iter + 1
        ! update
        if (mod(iter, 10) == 0) then
          call dscal(n, tmp, prevx1, 1)
          call update_omega(iter, n, prevx1, prevx2, omega)
        end if
      end do
      call daxpy(n, -1.0d0, x, incx, prevx, 1)

      rerror = maxval(abs(prevx / x(1,1:n)))

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovst_sor_csc

end module markovst_sor
