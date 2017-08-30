!
! mexp conv
!

module mexpconv_unif
  implicit none

contains

  ! Description: convolution integral operation for matrix exp form;
  !
  !                           |t
  ! transH(MH) = transH(MH) + | exp(transQ(Q)*s) * x * y' * exp(transQ(Q)*(t-s)) ds
  !                           |0
  !
  !        and
  !
  !        z = exp(transQ(Q)*t) * x
  !
  !        t is involved in the Poisson probability vector.
  !        qv is an uniformed parameter
  !        return value is z

  subroutine mexpconv_unif_dense_vec(transQ, transH, n, P, ldp, &
    qv, left, right, poi, weight, x, incx, y, incy, z, incz, H, ldh, atol)
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, ldp, left, right, incx, incy, incz, ldh
    double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), H(1:ldh,1:n)

    character :: itransQ
    integer :: k, l
    double precision :: xi(1:n), tmp(1:n)
    double precision, allocatable :: vc(:,:)
    allocate(vc(1:n, left:right))

    select case (transQ)
      case ('N', 'n')
        itransQ = 'T'
      case ('T', 't')
        itransQ = 'N'
    end select

    vc(1:n,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,right), 1)
    do l = right-1, left+1, -1
      call dgemv(itransQ, n, n, 1.0d0, P, ldp, vc(1,l+1), 1, 0.0d0, vc(1,l), 1)
      call daxpy(n, poi(l), y, incy, vc(1,l), 1)
    end do

    tmp(1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call dcopy(n, tmp, 1, z, incz)
    do k = 1, n
      call dscal(n, qv*weight, H(1,k), 1)
    end do
    call daxpy(n, poi(left), xi, 1, z, incz)

    select case (transH)
      case ('N', 'n')
        call dger(n, n, 1.0d0, xi, 1, vc(1,left+1), 1, H, ldh)
      case ('T', 't')
        call dger(n, n, 1.0d0, vc(1,left+1), 1, xi, 1, H, ldh)
    end select

    do l = left+1, right-1
      call dcopy(n, xi, 1, tmp, 1)
      call dgemv(transQ, n, n, 1.0d0, P, ldp, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(l), xi, 1, z, incz)
      select case (transH)
        case ('N', 'n')
          call dger(n, n, 1.0d0, xi, 1, vc(1,l+1), 1, H, ldh)
        case ('T', 't')
          call dger(n, n, 1.0d0, vc(1,l+1), 1, xi, 1, H, ldh)
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    do k = 1, n
      call dscal(n, 1.0/qv/weight, H(1,k), 1)
    end do

    deallocate(vc)
  end subroutine mexpconv_unif_dense_vec

  subroutine mexpconv_unif_csr_vec(transQ, transH, n, spP, rowptr, colind, nnz, &
    qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
    use spblas
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, nnz, left, right, incx, incy, incz
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

    character :: itransQ
    integer :: l
    double precision :: xi(1:n), tmp(1:n)
    double precision, allocatable :: vc(:,:)
    allocate(vc(1:n, left:right))

    select case (transQ)
      case ('N', 'n')
        itransQ = 'T'
      case ('T', 't')
        itransQ = 'N'
    end select

    vc(1:n,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,right), 1)
    do l = right-1, left+1, -1
      call dcsrmv(itransQ, n, n, 1.0d0, spP, rowptr, colind, nnz, vc(1,l+1), 1, 0.0d0, vc(1,l), 1)
      call daxpy(n, poi(l), y, incy, vc(1,l), 1)
    end do

    tmp(1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call dcopy(n, tmp, 1, z, incz)
    call dscal(nnz, qv*weight, spH, 1)
    call daxpy(n, poi(left), xi, 1, z, incz)

    select case (transH)
      case ('N', 'n')
        call dcsrr(n, n, 1.0d0, xi, 1, vc(1,left+1), 1, spH, rowptr, colind, nnz)
      case ('T', 't')
        call dcsrr(n, n, 1.0d0, vc(1,left+1), 1, xi, 1, spH, rowptr, colind, nnz)
    end select

    do l = left+1, right-1
      call dcopy(n, xi, 1, tmp, 1)
      call dcsrmv(transQ, n, n, 1.0d0, spP, rowptr, colind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(l), xi, 1, z, incz)
      select case (transH)
        case ('N', 'n')
          call dcsrr(n, n, 1.0d0, xi, 1, vc(1,l+1), 1, spH, rowptr, colind, nnz)
        case ('T', 't')
          call dcsrr(n, n, 1.0d0, vc(1,l+1), 1, xi, 1, spH, rowptr, colind, nnz)
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz, 1.0/qv/weight, spH, 1)

    deallocate(vc)
  end subroutine mexpconv_unif_csr_vec

  subroutine mexpconv_unif_csc_vec(transQ, transH, n, spP, colptr, rowind, nnz, &
    qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
    use spblas
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, nnz, left, right, incx, incy, incz
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

    character :: itransQ
    integer :: l
    double precision :: xi(1:n), tmp(1:n)
    double precision, allocatable :: vc(:,:)
    allocate(vc(1:n, left:right))

    select case (transQ)
      case ('N', 'n')
        itransQ = 'T'
      case ('T', 't')
        itransQ = 'N'
    end select

    vc(1:n,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,right), 1)
    do l = right-1, left+1, -1
      call dcscmv(itransQ, n, n, 1.0d0, spP, colptr, rowind, nnz, vc(1,l+1), 1, 0.0d0, vc(1,l), 1)
      call daxpy(n, poi(l), y, incy, vc(1,l), 1)
    end do

    tmp(1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call dcopy(n, tmp, 1, z, incz)
    call dscal(nnz, qv*weight, spH, 1)
    call daxpy(n, poi(left), xi, 1, z, incz)

    select case (transH)
      case ('N', 'n')
        call dcscr(n, n, 1.0d0, xi, 1, vc(1,left+1), 1, spH, colptr, rowind, nnz)
      case ('T', 't')
        call dcscr(n, n, 1.0d0, vc(1,left+1), 1, xi, 1, spH, colptr, rowind, nnz)
    end select

    do l = left+1, right-1
      call dcopy(n, xi, 1, tmp, 1)
      call dcscmv(transQ, n, n, 1.0d0, spP, colptr, rowind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(l), xi, 1, z, incz)
      select case (transH)
        case ('N', 'n')
          call dcscr(n, n, 1.0d0, xi, 1, vc(1,l+1), 1, spH, colptr, rowind, nnz)
        case ('T', 't')
          call dcscr(n, n, 1.0d0, vc(1,l+1), 1, xi, 1, spH, colptr, rowind, nnz)
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz, 1.0/qv/weight, spH, 1)

    deallocate(vc)
  end subroutine mexpconv_unif_csc_vec

  subroutine mexpconv_unif_coo_vec(transQ, transH, n, spP, rowind, colind, nnz, &
    qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
    use spblas
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, nnz, left, right, incx, incy, incz
    double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
    integer, intent(in) :: rowind(1:*), colind(1:*)
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

    character :: itransQ
    integer :: l
    double precision :: xi(1:n), tmp(1:n)
    double precision, allocatable :: vc(:,:)
    allocate(vc(1:n, left:right))

    select case (transQ)
      case ('N', 'n')
        itransQ = 'T'
      case ('T', 't')
        itransQ = 'N'
    end select

    vc(1:n,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,right), 1)
    do l = right-1, left+1, -1
      call dcoomv(itransQ, n, n, 1.0d0, spP, rowind, colind, nnz, vc(1,l+1), 1, 0.0d0, vc(1,l), 1)
      call daxpy(n, poi(l), y, incy, vc(1,l), 1)
    end do

    tmp(1:n) = 0.0d0
    call dcopy(n, x, incx, xi, 1)
    call dcopy(n, tmp, 1, z, incz)
    call dscal(nnz, qv*weight, spH, 1)
    call daxpy(n, poi(left), xi, 1, z, incz)

    select case (transH)
      case ('N', 'n')
        call dcoor(n, n, 1.0d0, xi, 1, vc(1,left+1), 1, spH, rowind, colind, nnz)
      case ('T', 't')
        call dcoor(n, n, 1.0d0, vc(1,left+1), 1, xi, 1, spH, rowind, colind, nnz)
    end select

    do l = left+1, right-1
      call dcopy(n, xi, 1, tmp, 1)
      call dcoomv(transQ, n, n, 1.0d0, spP, rowind, colind, nnz, tmp, 1, 0.0d0, xi, 1)
      call daxpy(n, poi(l), xi, 1, z, incz)
      select case (transH)
        case ('N', 'n')
          call dcoor(n, n, 1.0d0, xi, 1, vc(1,l+1), 1, spH, rowind, colind, nnz)
        case ('T', 't')
          call dcoor(n, n, 1.0d0, vc(1,l+1), 1, xi, 1, spH, rowind, colind, nnz)
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz, 1.0/qv/weight, spH, 1)

    deallocate(vc)
  end subroutine mexpconv_unif_coo_vec

end module mexpconv_unif

