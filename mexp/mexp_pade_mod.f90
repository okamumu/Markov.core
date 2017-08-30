!
! mexp pade
!

module mexp_pade
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

  subroutine mexp_pade_dense(trans, n, alpha, MA, lda, ME, lde, eps)
    character, intent(in) :: trans
    integer, intent(in) :: n, lda, lde
    double precision, intent(in) :: alpha, MA(1:lda,1:n), eps
    double precision, intent(out) :: ME(1:lde,1:n)

    double precision :: norma, tolerr, c
    integer :: i, j, k, q, info
    integer :: ipiv(1:n)
    double precision :: MN(1:n,1:n), MD(1:n,1:n), MX(1:n,1:n), MT(1:n,1:n)

    ! scaled A
    select case (trans)
      case ('N', 'n')
        ME(1:n,1:n) = alpha * MA(1:n,1:n)
      case ('T', 't')
        ME(1:n,1:n) = transpose(alpha * MA(1:n,1:n))
    end select
    norma = maxval(abs(ME(1:n,1:n)))
    j = max(0, 1 + floor(log(norma)/log(2.0d0)))
    ME(1:n,1:n) = ME(1:n,1:n) / 2**j

    ! find q
    q = 1
    tolerr = 1.0d0 / 6.0d0
    do while (tolerr > eps/norma)
      tolerr = tolerr / (16d0 * (3d0 + 4d0 * q * (2d0 + q)))
      q = q + 1
    end do

    ! initialize
    MD(1:n,1:n) = 0.0d0
    MN(1:n,1:n) = 0.0d0
    MX(1:n,1:n) = 0.0d0
    MT(1:n,1:n) = 0.0d0
    c = 1.0d0
    do i = 1,n
       MD(i,i) = 1.0d0
       MN(i,i) = 1.0d0
       MX(i,i) = 1.0d0
    end do

    ! make
    do k = 1, q
       c = c * real(q - k + 1, kind=8)/ real((2*q - k + 1) * k, kind=8)
       call dgemm('N', 'N', n, n, n, 1.0d0, ME, lde, MX, n, 0.0d0, MT, n)
       MX(1:n,1:n) = MT(1:n,1:n)
       call daxpy(n*n, c, MX, 1, MN, 1)
       if (mod(k,2) == 0) then
          call daxpy(n*n, c, MX, 1, MD, 1)
       else
          call daxpy(n*n, -c, MX, 1, MD, 1)
       end if
    end do
    ME(1:n,1:n) = MN(1:n,1:n)
    call dgesv(n, n, MD, n, ipiv, ME, lde, info)
    do k = 1, j
       call dgemm('N', 'N', n, n, n, 1.0d0, ME, lde, ME, lde, 0.0d0, MT, n)
       ME(1:n,1:n) = MT(1:n,1:n)
    end do
  end subroutine mexp_pade_dense

  subroutine mexp_pade_csr(trans, n, alpha, spMA, rowptr, colind, nnz, ME, lde, eps)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    double precision, intent(in) :: alpha, spMA(1:*), eps
    double precision, intent(out) :: ME(1:lde,1:n)

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call csr_to_dense(n, n, spMA, rowptr, colind, nnz, MA, n)
    call mexp_pade_dense(trans, n, alpha, MA, n, ME, lde, eps)

    deallocate(MA)

  end subroutine mexp_pade_csr

  subroutine mexp_pade_csc(trans, n, alpha, spMA, colptr, rowind, nnz, ME, lde, eps)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    double precision, intent(in) :: alpha, spMA(1:*), eps
    double precision, intent(out) :: ME(1:lde,1:n)

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call csc_to_dense(n, n, spMA, colptr, rowind, nnz, MA, n)
    call mexp_pade_dense(trans, n, alpha, MA, n, ME, lde, eps)

    deallocate(MA)

  end subroutine mexp_pade_csc

  subroutine mexp_pade_coo(trans, n, alpha, spMA, rowind, colind, nnz, ME, lde, eps)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde
    integer, intent(in) :: rowind(1:*), colind(1:*)
    double precision, intent(in) :: alpha, spMA(1:*), eps
    double precision, intent(out) :: ME(1:lde,1:n)

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call coo_to_dense(n, n, spMA, rowind, colind, nnz, MA, n)
    call mexp_pade_dense(trans, n, alpha, MA, n, ME, lde, eps)

    deallocate(MA)

  end subroutine mexp_pade_coo
end module mexp_pade

