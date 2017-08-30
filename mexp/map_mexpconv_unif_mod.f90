
module map_mexpconv_unif
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

  ! Description: convolution integral operation for matrix exp form;
  !
  !               u   |t
  !        MH0 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i)*(t-s)) ds
  !              i=0  |0
  !
  !              u-1  |t
  !        MH1 = Sum  | exp(D(i)*s) * vf * vb * exp(D(u-i-1)*(t-s)) ds
  !              i=0  |0
  !
  !        and
  !
  !        vf = vf * exp(D(u)*t)
  !
  !        where D is a MAP kernel, i.e,
  !
  !                 | Q0 Q1          |
  !                 |    Q0 Q1       |
  !          D(i) = |       Q0 Q1    |
  !                 |          .. .. |
  !                 |             Q0 |
  !
  !        The size of D(i) is given by i-by-i structured matrix,
  !        and t is involved in the Poisson probability vector.
  !
  ! Parameters
  !      n: size of CTMC (DTMC) kernel
  !      nnz: the number of non-zeros in CTMC kernel
  !      Q0: CTMC kernel (sparse matrix, phase matrix of MAP)
  !      rowptr0: row pointer vector of Q0
  !      colind0: column index vector of Q0
  !      Q1: CTMC kernel (sparse matrix, rate matrix of MAP)
  !      rowptr1: row pointer vector of Q1
  !      colind1: column index vector of Q1
  !      poi: Poisson probability
  !      right: right bound of Poisson range
  !      weight: normalizing constnt of a vector poi
  !      u: the number of arrivals
  !      vf: forward vector (inout)
  !      vb: backward vector (in)
  !      MH0: convint result (out), (sojourn time and phase transition)
  !      MH1: convint result (out), (arrial transition)

  subroutine map_mexpconv_unif_dense_transQ(transH, n, P0, ldp0, P1, ldp1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, ldh0, H1, ldh1, atol, xi, vc)
    character, intent(in) :: transH
    integer, intent(in) :: n, u, ldp0, ldp1, left, right, incx, incy, incz, ldh0, ldh1
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), H0(1:ldh0,1:n), H1(1:ldh1,1:n)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: k, j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,u,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dgemv('N', n, n, 1.0d0, P0, ldp0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= u) then
          call dgemv('N', n, n, 1.0d0, P1, ldp1, vc(1,j+1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,u,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,0), 1)
    do k = 1, n
      call dscal(n, qv*weight, H0(1,k), 1)
      call dscal(n, qv*weight, H1(1,k), 1)
    end do
    call daxpy(n, poi(left), xi(1,u), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, ldh0)
          if (j /= u) then
            call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,left+1), 1, H1, ldh1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dger(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, ldh0)
          if (j /= u) then
            call dger(n, n, 1.0d0, vc(1,j+1,left+1), 1, xi(1,j), 1, H1, ldh1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = u, 0, -1
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dgemv('T', n, n, 1.0d0, P0, ldp0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= 0) then
          call dgemv('T', n, n, 1.0d0, P1, ldp1, xi(1,j-1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,u), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, ldh0)
            if (j /= u) then
              call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,l+1), 1, H1, ldh1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dger(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, ldh0)
            if (j /= u) then
              call dger(n, n, 1.0d0, vc(1,j+1,l+1), 1, xi(1,j), 1, H1, ldh1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    do k = 1, n
      call dscal(n, 1.0/qv/weight, H0(1,k), 1)
      call dscal(n, 1.0/qv/weight, H1(1,k), 1)
    end do
  end subroutine map_mexpconv_unif_dense_transQ

  subroutine map_mexpconv_unif_dense_notransQ(transH, n, P0, ldp0, P1, ldp1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, ldh0, H1, ldh1, atol, xi, vc)
    character, intent(in) :: transH
    integer, intent(in) :: n, u, ldp0, ldp1, left, right, incx, incy, incz, ldh0, ldh1
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), H0(1:ldh0,1:n), H1(1:ldh1,1:n)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: k, j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,0,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dgemv('T', n, n, 1.0d0, P0, ldp0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= 0) then
          call dgemv('T', n, n, 1.0d0, P1, ldp1, vc(1,j-1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,0,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,u), 1)
    do k = 1, n
      call dscal(n, qv*weight, H0(1,k), 1)
      call dscal(n, qv*weight, H1(1,k), 1)
    end do
    call daxpy(n, poi(left), xi(1,0), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, ldh0)
          if (j /= 0) then
            call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,left+1), 1, H1, ldh1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dger(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, ldh0)
          if (j /= 0) then
            call dger(n, n, 1.0d0, vc(1,j-1,left+1), 1, xi(1,j), 1, H1, ldh1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = 0, u
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dgemv('N', n, n, 1.0d0, P0, ldp0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= u) then
          call dgemv('N', n, n, 1.0d0, P1, ldp1, xi(1,j+1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,0), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, ldh0)
            if (j /= 0) then
              call dger(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,l+1), 1, H1, ldh1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dger(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, ldh0)
            if (j /= 0) then
              call dger(n, n, 1.0d0, vc(1,j-1,l+1), 1, xi(1,j), 1, H1, ldh1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    do k = 1, n
      call dscal(n, 1.0/qv/weight, H0(1,k), 1)
      call dscal(n, 1.0/qv/weight, H1(1,k), 1)
    end do
  end subroutine map_mexpconv_unif_dense_notransQ

  subroutine map_mexpconv_unif_dense(transQ, transH, n, P0, ldp0, P1, ldp1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, ldh0, H1, ldh1, atol, xi, vc)
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, u, ldp0, ldp1, left, right, incx, incy, incz, ldh0, ldh1
    double precision, intent(in) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n), H0(1:ldh0,1:n), H1(1:ldh1,1:n)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    select case (transQ)
      case ('N', 'n')
        call map_mexpconv_unif_dense_notransQ(transH, n, P0, ldp0, P1, ldp1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, ldh0, H1, ldh1, atol, xi, vc)
      case ('T', 't')
        call map_mexpconv_unif_dense_transQ(transH, n, P0, ldp0, P1, ldp1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, ldh0, H1, ldh1, atol, xi, vc)
    end select
  end subroutine map_mexpconv_unif_dense

!!!!!!!!!!!!!!

  subroutine map_mexpconv_unif_csr_transQ(transH, n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    use spblas
    character, intent(in) :: transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,u,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dcsrmv('N', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= u) then
          call dcsrmv('N', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, vc(1,j+1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,u,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,0), 1)
    call dscal(nnz0, qv*weight, H0, 1)
    call dscal(nnz1, qv*weight, H1, 1)
    call daxpy(n, poi(left), xi(1,u), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, rowptr0, colind0, nnz0)
          if (j /= u) then
            call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,left+1), 1, H1, rowptr1, colind1, nnz1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dcsrr(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, rowptr0, colind0, nnz0)
          if (j /= u) then
            call dcsrr(n, n, 1.0d0, vc(1,j+1,left+1), 1, xi(1,j), 1, H1, rowptr1, colind1, nnz1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = u, 0, -1
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dcsrmv('T', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= 0) then
          call dcsrmv('T', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, xi(1,j-1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,u), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, rowptr0, colind0, nnz0)
            if (j /= u) then
              call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,l+1), 1, H1, rowptr1, colind1, nnz1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dcsrr(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, rowptr0, colind0, nnz0)
            if (j /= u) then
              call dcsrr(n, n, 1.0d0, vc(1,j+1,l+1), 1, xi(1,j), 1, H1, rowptr1, colind1, nnz1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz0, 1.0/qv/weight, H0, 1)
    call dscal(nnz1, 1.0/qv/weight, H1, 1)
  end subroutine map_mexpconv_unif_csr_transQ

  subroutine map_mexpconv_unif_csr_notransQ(transH, n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    use spblas
    character, intent(in) :: transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,0,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dcsrmv('T', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= 0) then
          call dcsrmv('T', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, vc(1,j-1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,0,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,u), 1)
    call dscal(nnz0, qv*weight, H0, 1)
    call dscal(nnz1, qv*weight, H1, 1)
    call daxpy(n, poi(left), xi(1,0), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, rowptr0, colind0, nnz0)
          if (j /= 0) then
            call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,left+1), 1, H1, rowptr1, colind1, nnz1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dcsrr(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, rowptr0, colind0, nnz0)
          if (j /= 0) then
            call dcsrr(n, n, 1.0d0, vc(1,j-1,left+1), 1, xi(1,j), 1, H1, rowptr1, colind1, nnz1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = 0, u
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dcsrmv('N', n, n, 1.0d0, P0, rowptr0, colind0, nnz0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= u) then
          call dcsrmv('N', n, n, 1.0d0, P1, rowptr1, colind1, nnz1, xi(1,j+1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,0), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, rowptr0, colind0, nnz0)
            if (j /= 0) then
              call dcsrr(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,l+1), 1, H1, rowptr1, colind1, nnz1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dcsrr(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, rowptr0, colind0, nnz0)
            if (j /= 0) then
              call dcsrr(n, n, 1.0d0, vc(1,j-1,l+1), 1, xi(1,j), 1, H1, rowptr1, colind1, nnz1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz0, 1.0/qv/weight, H0, 1)
    call dscal(nnz1, 1.0/qv/weight, H1, 1)
  end subroutine map_mexpconv_unif_csr_notransQ

  subroutine map_mexpconv_unif_csr(transQ, transH, n, P0, rowptr0, colind0, nnz0, &
    P1, rowptr1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    select case (transQ)
      case ('N', 'n')
        call map_mexpconv_unif_csr_notransQ(transH, n, P0, rowptr0, colind0, nnz0, &
          P1, rowptr1, colind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
      case ('T', 't')
        call map_mexpconv_unif_csr_transQ(transH, n, P0, rowptr0, colind0, nnz0, &
          P1, rowptr1, colind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
    end select
  end subroutine map_mexpconv_unif_csr

  subroutine map_mexpconv_unif_csc(transQ, transH, n, P0, colptr0, rowind0, nnz0, &
    P1, colptr1, rowind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
    integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    character :: itransH
    select case (transH)
      case ('N', 'n')
        itransH = 'T'
      case ('T', 't')
        itransH = 'N'
    end select

    select case (transQ)
      case ('T', 't')
        call map_mexpconv_unif_csr_notransQ(itransH, n, P0, colptr0, rowind0, nnz0, &
          P1, colptr1, rowind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
      case ('N', 'n')
        call map_mexpconv_unif_csr_transQ(itransH, n, P0, colptr0, rowind0, nnz0, &
          P1, colptr1, rowind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
    end select
  end subroutine map_mexpconv_unif_csc

  subroutine map_mexpconv_unif_coo_transQ(transH, n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    use spblas
    character, intent(in) :: transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,u,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dcoomv('N', n, n, 1.0d0, P0, rowind0, colind0, nnz0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= u) then
          call dcoomv('N', n, n, 1.0d0, P1, rowind1, colind1, nnz1, vc(1,j+1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,u,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,0), 1)
    call dscal(nnz0, qv*weight, H0, 1)
    call dscal(nnz1, qv*weight, H1, 1)
    call daxpy(n, poi(left), xi(1,u), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, rowind0, colind0, nnz0)
          if (j /= u) then
            call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,left+1), 1, H1, rowind1, colind1, nnz1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dcoor(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, rowind0, colind0, nnz0)
          if (j /= u) then
            call dcoor(n, n, 1.0d0, vc(1,j+1,left+1), 1, xi(1,j), 1, H1, rowind1, colind1, nnz1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = u, 0, -1
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dcoomv('T', n, n, 1.0d0, P0, rowind0, colind0, nnz0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= 0) then
          call dcoomv('T', n, n, 1.0d0, P1, rowind1, colind1, nnz1, xi(1,j-1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,u), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, rowind0, colind0, nnz0)
            if (j /= u) then
              call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j+1,l+1), 1, H1, rowind1, colind1, nnz1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dcoor(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, rowind0, colind0, nnz0)
            if (j /= u) then
              call dcoor(n, n, 1.0d0, vc(1,j+1,l+1), 1, xi(1,j), 1, H1, rowind1, colind1, nnz1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz0, 1.0/qv/weight, H0, 1)
    call dscal(nnz1, 1.0/qv/weight, H1, 1)
  end subroutine map_mexpconv_unif_coo_transQ

  subroutine map_mexpconv_unif_coo_notransQ(transH, n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    use spblas
    character, intent(in) :: transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    integer :: j, l
    double precision :: tmp(1:n)

    vc(1:n,0:u,right) = 0.0d0
    call daxpy(n, poi(right), y, incy, vc(1,0,right), 1)
    do l = right-1, left+1, -1
      do j = 0, u
        call dcoomv('T', n, n, 1.0d0, P0, rowind0, colind0, nnz0, vc(1,j,l+1), 1, 0.0d0, vc(1,j,l), 1)
        if (j /= 0) then
          call dcoomv('T', n, n, 1.0d0, P1, rowind1, colind1, nnz1, vc(1,j-1,l+1), 1, 1.0d0, vc(1,j,l), 1)
        end if
      end do
      call daxpy(n, poi(l), y, incy, vc(1,0,l), 1)
    end do

    z(1,1:n) = 0.0d0
    xi(1:n,0:u) = 0.0d0
    call dcopy(n, x, incx, xi(1,u), 1)
    call dscal(nnz0, qv*weight, H0, 1)
    call dscal(nnz1, qv*weight, H1, 1)
    call daxpy(n, poi(left), xi(1,0), 1, z, incz)

    select case (transH)
      case ('N', 'n')
        do j = 0, u
          call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j,left+1), 1, H0, rowind0, colind0, nnz0)
          if (j /= 0) then
            call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,left+1), 1, H1, rowind1, colind1, nnz1)
          end if
        end do
      case ('T', 't')
        do j = 0, u
          call dcoor(n, n, 1.0d0, vc(1,j,left+1), 1, xi(1,j), 1, H0, rowind0, colind0, nnz0)
          if (j /= 0) then
            call dcoor(n, n, 1.0d0, vc(1,j-1,left+1), 1, xi(1,j), 1, H1, rowind1, colind1, nnz1)
          end if
        end do
    end select

    do l = left+1, right-1
      do j = 0, u
        call dcopy(n, xi(1,j), 1, tmp, 1)
        call dcoomv('N', n, n, 1.0d0, P0, rowind0, colind0, nnz0, tmp, 1, 0.0d0, xi(1,j), 1)
        if (j /= u) then
          call dcoomv('N', n, n, 1.0d0, P1, rowind1, colind1, nnz1, xi(1,j+1), 1, 1.0d0, xi(1,j), 1)
        end if
      end do
      call daxpy(n, poi(l), xi(1,0), 1, z, incz)
      select case (transH)
        case ('N', 'n')
          do j = 0, u
            call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j,l+1), 1, H0, rowind0, colind0, nnz0)
            if (j /= 0) then
              call dcoor(n, n, 1.0d0, xi(1,j), 1, vc(1,j-1,l+1), 1, H1, rowind1, colind1, nnz1)
            end if
          end do
        case ('T', 't')
          do j = 0, u
            call dcoor(n, n, 1.0d0, vc(1,j,l+1), 1, xi(1,j), 1, H0, rowind0, colind0, nnz0)
            if (j /= 0) then
              call dcoor(n, n, 1.0d0, vc(1,j-1,l+1), 1, xi(1,j), 1, H1, rowind1, colind1, nnz1)
            end if
          end do
      end select
    end do
    call dscal(n, 1.0/weight, z, incz)
    call dscal(nnz0, 1.0/qv/weight, H0, 1)
    call dscal(nnz1, 1.0/qv/weight, H1, 1)
  end subroutine map_mexpconv_unif_coo_notransQ

  subroutine map_mexpconv_unif_coo(transQ, transH, n, P0, rowind0, colind0, nnz0, &
    P1, rowind1, colind1, nnz1, &
    qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
    H0, H1, atol, xi, vc)
    character, intent(in) :: transQ, transH
    integer, intent(in) :: n, u, nnz0, nnz1, left, right, incx, incy, incz
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: qv, weight, poi(left:right), atol
    double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
    double precision, intent(out) :: z(1:incz,1:n)
    double precision, intent(out) :: H0(1:nnz0), H1(1:nnz1)
    double precision :: xi(1:n,0:u)
    double precision :: vc(1:n,0:u,left:right)

    select case (transQ)
      case ('N', 'n')
        call map_mexpconv_unif_coo_notransQ(transH, n, P0, rowind0, colind0, nnz0, &
          P1, rowind1, colind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
      case ('T', 't')
        call map_mexpconv_unif_coo_transQ(transH, n, P0, rowind0, colind0, nnz0, &
          P1, rowind1, colind1, nnz1, &
          qv, left, right, poi, weight, u, x, incx, y, incy, z, incz, &
          H0, H1, atol, xi, vc)
    end select
  end subroutine map_mexpconv_unif_coo

end module map_mexpconv_unif

