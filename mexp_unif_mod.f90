!
! mexp pade
!

module mexp_unif
  implicit none

contains

  !
  ! Description:
  !
  !   ME = exp(alpha trans(MA))
  !
  ! Parameters:
  !   n: size of square matrix A
  !

  subroutine mexp_unif_dense_vec(trans, n, P, ldp, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    character, intent(in) :: trans
    integer, intent(in) :: n, ldp, left, right, incx, incy
    double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), x(1:*), atol
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: xi(1:n), tmp(1:n)

    call dcopy(n, x, incx, xi, 1)
    y(1,1:n) = 0.0d0
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcopy(n, xi, 1, tmp, 1)
      call dgemv(trans, n, n, 1.0d0, P, ldp, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
      if (sum(xi) < atol) then
        exit
      end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine mexp_unif_dense_vec

  subroutine mexp_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, incx, incy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: xi(1:n), tmp(1:n)

    call dcopy(n, x, incx, xi, 1)
    y(1,1:n) = 0.0d0
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcopy(n, xi, 1, tmp, 1)
      call dcsrmv(trans, n, n, 1.0d0, spP, rowptr, colind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
      if (sum(xi) < atol) then
        exit
      end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine mexp_unif_csr_vec

  subroutine mexp_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, incx, incy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: xi(1:n), tmp(1:n)

    call dcopy(n, x, incx, xi, 1)
    y(1,1:n) = 0.0d0
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcopy(n, xi, 1, tmp, 1)
      call dcscmv(trans, n, n, 1.0d0, spP, colptr, rowind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
      if (sum(xi) < atol) then
        exit
      end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine mexp_unif_csc_vec

  subroutine mexp_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
    left, right, poi, weight, x, incx, y, incy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, incx, incy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
    integer, intent(in) :: rowind(1:*), colind(1:*)
    double precision, intent(out) :: y(1:incy,1:n)

    integer :: k
    double precision :: xi(1:n), tmp(1:n)

    call dcopy(n, x, incx, xi, 1)
    y(1,1:n) = 0.0d0
    call daxpy(n, poi(left), xi, 1, y, incy)
    do k = left+1, right
      call dcopy(n, xi, 1, tmp, 1)
      call dcoomv(trans, n, n, 1.0d0, spP, rowind, colind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(k), xi, 1, y, incy)
      if (sum(xi) < atol) then
        exit
      end if
    end do
    call dscal(n, 1.0d0/weight, y, incy)
  end subroutine mexp_unif_coo_vec

  !
  ! Description:
  !
  !   ME = exp(alpha trans(MA))
  !
  ! Parameters:
  !   n: size of square matrix A
  !

  subroutine mexp_unif_dense_mat(trans, n, P, ldp, qv, &
    left, right, poi, weight, m, x, ldx, y, ldy, atol)
    character, intent(in) :: trans
    integer, intent(in) :: n, ldp, left, right, m, ldx, ldy
    double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:ldx,1:m)
    double precision, intent(out) :: y(1:ldy,1:m)

    integer :: j, k
    double precision :: xi(1:n,1:m), tmp(1:n,1:m)

    xi(1:n,1:m) = x(1:n,1:m)
    y(1:n,1:m) = 0.0d0
    do j = 1, m
      call daxpy(n, poi(left), xi(1,j), 1, y(1,j), 1)
    end do
    do k = left+1, right
      call dcopy(n*m, xi, 1, tmp, 1)
      call dgemm(trans, 'N', n, m, n, 1.0d0, P, ldp, tmp, n, 0.0d0, xi, n)
      do j = 1, m
        call daxpy(n, poi(k), xi(1,j), 1, y(1,j), 1)
      end do
      if (sum(xi) < atol) then
        exit
      end if
    end do
    do j = 1, m
      call dscal(n, 1.0d0/weight, y(1,j), 1)
    end do
  end subroutine mexp_unif_dense_mat

  subroutine mexp_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
    left, right, poi, weight, m, x, ldx, y, ldy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    double precision, intent(in) :: x(1:ldx,1:m)
    double precision, intent(out) :: y(1:ldy,1:m)

    integer :: j, k
    double precision :: xi(1:n,1:m), tmp(1:n,1:m)

    xi(1:n,1:m) = x(1:n,1:m)
    y(1:n,1:m) = 0.0d0
    do j = 1, m
      call daxpy(n, poi(left), xi(1,j), 1, y(1,j), 1)
    end do
    do k = left+1, right
      call dcopy(n*m, xi, 1, tmp, 1)
      call dcsrmm(trans, 'N', n, m, n, 1.0d0, spP, rowptr, colind, nnz, tmp, n, 0.0d0, xi, n)
      do j = 1, m
        call daxpy(n, poi(k), xi(1,j), 1, y(1,j), 1)
      end do
      if (sum(xi) < atol) then
        exit
      end if
    end do
    do j = 1, m
      call dscal(n, 1.0d0/weight, y(1,j), 1)
    end do
  end subroutine mexp_unif_csr_mat

  subroutine mexp_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
    left, right, poi, weight, m, x, ldx, y, ldy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    double precision, intent(in) :: x(1:ldx,1:m)
    double precision, intent(out) :: y(1:ldy,1:m)

    integer :: j, k
    double precision :: xi(1:n,1:m), tmp(1:n,1:m)

    xi(1:n,1:m) = x(1:n,1:m)
    y(1:n,1:m) = 0.0d0
    do j = 1, m
      call daxpy(n, poi(left), xi(1,j), 1, y(1,j), 1)
    end do
    do k = left+1, right
      call dcopy(n*m, xi, 1, tmp, 1)
      call dcscmm(trans, 'N', n, m, n, 1.0d0, spP, colptr, rowind, nnz, tmp, n, 0.0d0, xi, n)
      do j = 1, m
        call daxpy(n, poi(k), xi(1,j), 1, y(1,j), 1)
      end do
      if (sum(xi) < atol) then
        exit
      end if
    end do
    do j = 1, m
      call dscal(n, 1.0d0/weight, y(1,j), 1)
    end do
  end subroutine mexp_unif_csc_mat

  subroutine mexp_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
    left, right, poi, weight, m, x, ldx, y, ldy, atol)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: rowind(1:*), colind(1:*)
    double precision, intent(in) :: x(1:ldx,1:m)
    double precision, intent(out) :: y(1:ldy,1:m)

    integer :: j, k
    double precision :: xi(1:n,1:m), tmp(1:n,1:m)

    xi(1:n,1:m) = x(1:n,1:m)
    y(1:n,1:m) = 0.0d0
    do j = 1, m
      call daxpy(n, poi(left), xi(1,j), 1, y(1,j), 1)
    end do
    do k = left+1, right
      call dcopy(n*m, xi, 1, tmp, 1)
      call dcoomm(trans, 'N', n, m, n, 1.0d0, spP, rowind, colind, nnz, tmp, n, 0.0d0, xi, n)
      do j = 1, m
        call daxpy(n, poi(k), xi(1,j), 1, y(1,j), 1)
      end do
      if (sum(xi) < atol) then
        exit
      end if
    end do
    do j = 1, m
      call dscal(n, 1.0d0/weight, y(1,j), 1)
    end do
  end subroutine mexp_unif_coo_mat

end module mexp_unif

