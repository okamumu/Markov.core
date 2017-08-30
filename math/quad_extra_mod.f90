!
! accerelation
!

module accerelation
  implicit none
  public aitken

contains

! Aitken's delta-squared accerelation
!
!   y = x0 - (Delta x0)^2 / Delta^2 x0
!     = x2 - (x2 - x1)^2 / (x2 - 2 x1 + x0)
!     = x0 - (x1 - x0)^2 / (x2 - 2 x1 + x0)
!
!         x0: a vector at the n-th step
!         x1: a vector at the n+1-st step
!         x2: a vector at the n+2-nd step
!         y: accelerated vector
!

  ! subroutine aitken(n, x0, incx0, x1, incx1, x2, incx2, y, incy)
  !   integer, intent(in) :: n, incx0, incx1, incx2, incy
  !   double precision, intent(in) :: x0(1:incx0,1:n), x1(1:incx1,1:n), x2(1:incx2,1:n)
  !   double precision, intent(out) :: y(1:incy,1:n)
  !
  !   ! double precision :: g(1:n), h(1:n)
  !   !
  !   ! g(1:n) = x0(1,1:n)
  !   ! daxpy(n, -1.0d0, x1, incx1, g, 1)
  !   !
  !   ! h(1:n) = x2(1,1:n)
  !   ! daxpy(n, -2.0d0, x1, incx1, h, 1)
  !   ! daxpy(n, 1.0d0, x0, incx0, h, 1)
  !   !
  !   ! g(1:n) = g(1:n) * g(1:n) / h(1:n)
  !   ! y(1,1:n) = x(1,1:n)
  !   ! daxpy(n, -1.0d0, g, 1, y, incy)
  !
  !   y(1,1:n) = x0(1,1:n) - (x1(1,1:n) - x0(1,1:n))^2 / (x2(1,1:n) - 2.0d0 * x1(1,1:n) + x0(1,1:n))
  ! end subroutine

  subroutine aitken(n, x0, incx0, x1, incx1, x2, incx2, y, incy)
    integer, intent(in) :: n, incx0, incx1, incx2, incy
    double precision, intent(in) :: x0(1:incx0,1:n), x1(1:incx1,1:n), x2(1:incx2,1:n)
    double precision, intent(out) :: y(1:incy,1:n)

    ! double precision :: g(1:n), h(1:n)
    !
    ! g(1:n) = x0(1,1:n)
    ! daxpy(n, -1.0d0, x1, incx1, g, 1)
    !
    ! h(1:n) = x2(1,1:n)
    ! daxpy(n, -2.0d0, x1, incx1, h, 1)
    ! daxpy(n, 1.0d0, x0, incx0, h, 1)
    !
    ! g(1:n) = g(1:n) * g(1:n) / h(1:n)
    ! y(1,1:n) = x(1,1:n)
    ! daxpy(n, -1.0d0, g, 1, y, incy)

    y(1,1:n) = x0(1,1:n) - (x1(1,1:n) - x0(1,1:n))^2 / (x2(1,1:n) - 2.0d0 * x1(1,1:n) + x0(1,1:n))
  end subroutine

end module accerelation
