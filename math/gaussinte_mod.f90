
module gaussinte
  implicit none
  double precision, private, parameter :: pai = 4.0d0 * atan(1.0d0)
  private
  public gaussinte_w, gaussinte_fx, gaussinte_fv

contains

!   Description: Gauss quadrature for the following integral

!      | b
!      |  f(x) dx
!      | a

!   gaussinte_w: make points and weights for n discrete points
!     n (in): the number of points. This is the size of both x and w.
!     x (out): x-axis points in the interval [-1, 1].
!     w (out): weights for the points.
!     eps (in): tolerance error.

!   gaussinte_fx: make x points for the interval [a,b]
!     n (in): the number of points.
!     x (in): the x points for interval [-1, 1].
!     a, b (in): lower and upper points for the integral
!     fx (out): x point for the interval [a, b]
!     return value: (b-a)/2

!   gaussinte_fv: compute the integral
!     n (in): the number of points.
!     w (in): weights for the x points in [-1,1]
!     c (in): (b-a)/2 ?
!     fv (in): function values at x points derived by gauss_inte_fx
!     return value: the interal value

  subroutine gaussinte_w(n, x, w, eps)
    integer, intent(in) :: n
    double precision, intent(out) :: x(1:n), w(1:n)
    double precision, intent(in) :: eps

    integer :: m, i, l
    double precision :: npai, tmp, dt
    double precision :: p0, p1, p2
    double precision :: q0, q1, q2

    select case (n)
      case (1)
        x(1) = 0.0d0
        w(1) = 2.0d0
        return
      case (2)
        x(1) = -sqrt(1.0d0/3.0d0)
        w(1) = 1.0d0
        x(2) = -x(1)
        w(2) = w(1)
        return
      case (3)
        x(1) = -sqrt(0.6d0)
        w(1) = 5.0d0/9.0d0
        x(2) = 0.0d0
        w(2) = 8.0d0/9.0d0
        x(3) = -x(1)
        w(3) = w(1)
        return
    end select

    m = n / 2
    npai = pai / (dble(n) + 0.5d0)
    do i = 1, m
      tmp = cos((i - 0.25d0) * npai)
      dt = tmp
      do while (abs(dt) > abs(tmp) * eps)
        p1 = tmp
        p2 = (3.0d0 * tmp**2 - 1.0d0) * 0.5d0
        q1 = 1.0d0
        q2 = 3.0d0 * tmp
        do l = 3, n
          p0 = p1
          p1 = p2
          p2 = ((l + l - 1) * tmp * p1 - (l-1) * p0) / l
          q0 = q1
          q1 = q2
          q2 = ((l + l - 1) * (tmp * q1 + p1) - (l-1) * q0) / l
        end do
        dt = p2 / q2
        tmp = tmp - dt
      end do
      x(i) = tmp
      w(i) = 2.0d0 / (n * p1 * q2)
    end do

    if (mod(n, 2) == 1) then
      tmp = dble(n)
      do i = 1, m
        tmp = tmp * (0.5d0 - i) / dble(i)
      end do
      x(m+1) = 0.0d0
      w(m+1) = 2.0d0 / tmp**2
    end if

    do i = 1, m
      x(n+1-i) = x(i)
      x(i) = -x(i) ! reverse order
      w(n+1-i) = w(i)
    end do
  end subroutine gaussinte_w

  function gaussinte_fx(n, x, a, b, fx) result(t1)
    integer, intent(in) :: n
    double precision, intent(in) :: x(1:n), a, b
    double precision, intent(out) :: fx(1:n)

    integer :: i
    double precision :: t1, t2

    t1 = (b - a)/2.0d0
    t2 = (b + a)/2.0d0
    do i = 1, n
      fx(i) = t1 * x(i) + t2
    end do
  end function gaussinte_fx

  function gaussinte_fv(n, w, c, fv) result(s)
    integer, intent(in) :: n
    double precision, intent(in) :: w(1:n), c, fv(1:n)

    integer :: i
    double precision :: s

    s = 0.0d0
    do i = 1, n
      s = s + w(i) * fv(i)
    end do
    s = s * c
  end function gaussinte_fv

end module gaussinte
