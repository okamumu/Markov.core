!
! mexp pade
!

module map_mexp_unif
  implicit none
  private map_mexp_unif_dense_notrans
  private map_mexp_unif_dense_trans
  private map_mexp_unif_dense_NA
  private map_mexp_unif_csr_notrans
  private map_mexp_unif_csr_trans
  private map_mexp_unif_csr_NA
  public map_mexp_unif_dense, map_mexp_unif_csr
contains

!     // ! Description: vector-matrix operation for matrix exp form;
!     // !
!     // !        y = exp(D*t) * x
!     // !
!     // !        where D is a MAP kernel, i.e,
!     // !
!     // !              | Q0 Q1          |
!     // !              |    Q0 Q1       |
!     // !          D = |       Q0 Q1    |
!     // !              |          .. .. |
!     // !              |             Q0 |
!     // !
!     // !        The size of D is given by u-by-u structured matrix,
!     // !        and t is involved in the Poisson probability vector.
!     // !
!     // ! Parameters
!     // !      n: size of CTMC (DTMC) kernel
!     // !      nnz: the number of non-zeros in CTMC kernel
!     // !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
!     // !      rowptr0: row pointer vector of Q0
!     // !      colind0: column index vector of Q0
!     // !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
!     // !      rowptr1: row pointer vector of Q1
!     // !      colind1: column index vector of Q1
!     // !      poi: Poisson probability
!     // !      right: right bound of Poisson range
!     // !      weight: normalizing constnt of a vector poi
!     // !      u: the number of arrivals
!     // !      x: input vector
!     // !      y: output vector

! private

  subroutine map_mexp_unif_dense_notrans(n, P0, ldp0, P1, ldp1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    integer, intent(in) :: n, ldp0, ldp1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n), x(1:incx,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,u) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,0), 1, y, incy)
    do k = left+1, right
      do l = 0, u
        call dgemv('N', n, n, 1.0d0, P0, ldp0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= u) then
          call dgemv('N', n, n, 1.0d0, P1, ldp1, xi(1,l+1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,0), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_dense_notrans

  subroutine map_mexp_unif_dense_trans(n, P0, ldp0, P1, ldp1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    integer, intent(in) :: n, ldp0, ldp1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n), x(1:incx,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,0) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,u), 1, y, incy)
    do k = left+1, right
      do l = u, 0, -1
        call dgemv('T', n, n, 1.0d0, P0, ldp0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= 0) then
          call dgemv('T', n, n, 1.0d0, P1, ldp1, xi(1,l-1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,u), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_dense_trans

  subroutine map_mexp_unif_dense_NA(trans, n, P0, ldp0, P1, ldp1, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    character, intent(in) :: trans
    integer, intent(in) :: n, ldp0, ldp1, left, right, incx, incy
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n), x(1:incx,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: tmp(1:n), xi(1:n)

    y(1,1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dgemv(trans, n, n, 1.0d0, P0, ldp0, xi, 1, 0.0d0, tmp, 1)
      call dgemv(trans, n, n, 1.0d0, P1, ldp1, xi, 1, 1.0d0, tmp, 1)
      call dcopy(n, tmp, 1, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_dense_NA

! public

  subroutine map_mexp_unif_dense(trans, n, P0, ldp0, P1, ldp1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    character, intent(in) :: trans
    integer, intent(in) :: n, ldp0, ldp1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n), x(1:incx,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    if (u >= 0) then
      select case (trans)
        case ('N', 'n')
          call map_mexp_unif_dense_notrans(n, P0, ldp0, P1, ldp1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
        case ('T', 't')
          call map_mexp_unif_dense_trans(n, P0, ldp0, P1, ldp1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
      end select
    else
      call map_mexp_unif_dense_NA(trans, n, P0, ldp0, P1, ldp1, qv, &
        left, right, poi, weight, x, incx, y, incy, atol)
    end if
  end subroutine map_mexp_unif_dense

!!!!!!!!!!!!!!!!!!

  subroutine map_mexp_unif_csr_notrans(n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,u) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,0), 1, y, incy)
    do k = left+1, right
      do l = 0, u
        call dcsrmv('N', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= u) then
          call dcsrmv('N', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, xi(1,l+1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,0), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_csr_notrans

  subroutine map_mexp_unif_csr_trans(n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,0) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,u), 1, y, incy)
    do k = left+1, right
      do l = u, 0, -1
        call dcsrmv('T', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= 0) then
          call dcsrmv('T', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, xi(1,l-1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,u), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_csr_trans

  subroutine map_mexp_unif_csr_NA(trans, n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz0, nnz1, left, right, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: tmp(1:n), xi(1:n)

    y(1,1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcsrmv(trans, n, n, 1.0d0, P0, rowptr0, colind0, nnz0, xi, 1, 0.0d0, tmp, 1)
      call dcsrmv(trans, n, n, 1.0d0, P1, rowptr1, colind1, nnz1, xi, 1, 1.0d0, tmp, 1)
      call dcopy(n, tmp, 1, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_csr_NA

  subroutine map_mexp_unif_csr(trans, n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    if (u >= 0) then
      select case (trans)
        case ('N', 'n')
          call map_mexp_unif_csr_notrans(n, P0, rowptr0, colind0, nnz0, &
            P1, rowptr1, colind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
        case ('T', 't')
          call map_mexp_unif_csr_trans(n, P0, rowptr0, colind0, nnz0, &
            P1, rowptr1, colind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
      end select
    else
      call map_mexp_unif_csr_NA(trans, n, P0, rowptr0, colind0, nnz0, &
        P1, rowptr1, colind1, nnz1, qv, &
        left, right, poi, weight, x, incx, y, incy, atol)
    end if
  end subroutine map_mexp_unif_csr

  subroutine map_mexp_unif_csc(trans, n, P0, colptr0, rowind0, nnz0, &
    P1, colptr1, rowind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
    integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    if (u >= 0) then
      select case (trans)
        case ('N', 'n')
          call map_mexp_unif_csr_trans(n, P0, colptr0, rowind0, nnz0, &
            P1, colptr1, rowind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
        case ('T', 't')
          call map_mexp_unif_csr_notrans(n, P0, colptr0, rowind0, nnz0, &
            P1, colptr1, rowind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
      end select
    else
      select case (trans)
        case ('N', 'n')
          call map_mexp_unif_csr_NA('T', n, P0, colptr0, rowind0, nnz0, &
            P1, colptr1, rowind1, nnz1, qv, &
            left, right, poi, weight, x, incx, y, incy, atol)
        case ('T', 't')
          call map_mexp_unif_csr_NA('N', n, P0, colptr0, rowind0, nnz0, &
            P1, colptr1, rowind1, nnz1, qv, &
            left, right, poi, weight, x, incx, y, incy, atol)
      end select
    end if
  end subroutine map_mexp_unif_csc

  subroutine map_mexp_unif_coo_notrans(n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowind0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,u) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,0), 1, y, incy)
    do k = left+1, right
      do l = 0, u
        call dcoomv('N', n, n, 1.0d0, P0, rowind0, colind0, nnz0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= u) then
          call dcoomv('N', n, n, 1.0d0, P1, rowind1, colind1, nnz1, xi(1,l+1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,0), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_coo_notrans

  subroutine map_mexp_unif_coo_trans(n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowind0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    double precision :: tmp(1:n)
    integer :: k, l

    y(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    xi(1:n,0) = x(1,1:n)

    call daxpy(n, poi(left), xi(1,u), 1, y, incy)
    do k = left+1, right
      do l = u, 0, -1
        call dcoomv('T', n, n, 1.0d0, P0, rowind0, colind0, nnz0, xi(1,l), 1, 0.0d0, tmp, 1)
        if (l /= 0) then
          call dcoomv('T', n, n, 1.0d0, P1, rowind1, colind1, nnz1, xi(1,l-1), 1, 1.0d0, tmp, 1)
        end if
        call dcopy(n, tmp, 1, xi(1,l), 1)
      end do
      call daxpy(n, poi(k), xi(1,u), 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_coo_trans

  subroutine map_mexp_unif_coo_NA(trans, n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz0, nnz1, left, right, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowind0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: tmp(1:n), xi(1:n)

    y(1,1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcoomv(trans, n, n, 1.0d0, P0, rowind0, colind0, nnz0, xi, 1, 0.0d0, tmp, 1)
      call dcoomv(trans, n, n, 1.0d0, P1, rowind1, colind1, nnz1, xi, 1, 1.0d0, tmp, 1)
      call dcopy(n, tmp, 1, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
!       if (sum(xi) < atol) then
!         exit
!       end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine map_mexp_unif_coo_NA

  subroutine map_mexp_unif_coo(trans, n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, qv, &
    left, right, poi, weight, u, x, incx, y, incy, atol, xi)
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz0, nnz1, left, right, u, incx, incy
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1), x(1:incx,1:n)
    integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(out) :: y(1:incy,1:n)
    double precision :: xi(1:n,0:u)

    if (u >= 0) then
      select case (trans)
        case ('N', 'n')
          call map_mexp_unif_coo_notrans(n, P0, rowind0, colind0, nnz0, &
            P1, rowind1, colind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
        case ('T', 't')
          call map_mexp_unif_coo_trans(n, P0, rowind0, colind0, nnz0, &
            P1, rowind1, colind1, nnz1, qv, &
            left, right, poi, weight, u, x, incx, y, incy, atol, xi)
      end select
    else
      call map_mexp_unif_coo_NA(trans, n, P0, rowind0, colind0, nnz0, &
        P1, rowind1, colind1, nnz1, qv, &
        left, right, poi, weight, x, incx, y, incy, atol)
    end if
  end subroutine map_mexp_unif_coo
end module map_mexp_unif

