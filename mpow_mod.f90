
module mpow
  implicit none

contains

  !
  ! Description:
  !
  !   ME = trans(MA)^m
  !
  !

  subroutine mpow_dense(trans, n, MA, lda, ME, lde, m, info)
    character, intent(in) :: trans
    integer, intent(in) :: n, lda, lde, m
    double precision, intent(in) :: MA(1:lda,1:n)
    double precision, intent(out) :: ME(1:lde,1:n)
    integer, intent(out) :: info

    integer :: i, y
    integer :: ipiv(1:n)
    double precision :: MX(1:n,1:n), MT(1:n,1:n)

    info = 0
    if (m < 0) then
      select case (trans)
        case ('N', 'n')
          MT(1:n,1:n) = MA(1:n,1:n)
        case ('T', 't')
          MT(1:n,1:n) = transpose(MA(1:n,1:n))
      end select

      MX(1:n,1:n) = 0.0d0
      do i = 1, n
        MX(i,i) = 1.0d0
      end do
      call dgesv(n, n, MT, n, ipiv, MX, n, info)
      y = -m
    else
      select case (trans)
        case ('N', 'n')
          MX(1:n,1:n) = MA(1:n,1:n)
        case ('T', 't')
          MX(1:n,1:n) = transpose(MA(1:n,1:n))
      end select
      y = m
    end if

    ME(1:n,1:n) = 0.0d0
    do i = 1, n
      ME(i,i) = 1.0d0
    end do

    do while (y /= 0)
      if (mod(y,2) == 1) then
        MT(1:n,1:n) = ME(1:n,1:n)
        call dgemm('N', 'N', n, n, n, 1.0d0, MX, n, MT, n, 0.0d0, ME, lde)
      end if
      y = y / 2
      call dgemm('N', 'N', n, n, n, 1.0d0, MX, n, MX, n, 0.0d0, MT, n)
      MX(1:n,1:n) = MT(1:n,1:n)
    end do
  end subroutine mpow_dense

  subroutine mpow_csr(trans, n, spMA, rowptr, colind, nnz, ME, lde, m, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde, m
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    double precision, intent(in) :: spMA(1:*)
    double precision, intent(out) :: ME(1:lde,1:n)
    integer, intent(out) :: info

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call csr_to_dense(n, n, spMA, rowptr, colind, nnz, MA, n)
    call mpow_dense(trans, n, MA, n, ME, lde, m, info)

    deallocate(MA)
  end subroutine mpow_csr

  subroutine mpow_csc(trans, n, spMA, colptr, rowind, nnz, ME, lde, m, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde, m
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    double precision, intent(in) :: spMA(1:*)
    double precision, intent(out) :: ME(1:lde,1:n)
    integer, intent(out) :: info

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call csc_to_dense(n, n, spMA, colptr, rowind, nnz, MA, n)
    call mpow_dense(trans, n, MA, n, ME, lde, m, info)

    deallocate(MA)
  end subroutine mpow_csc

  subroutine mpow_coo(trans, n, spMA, rowind, colind, nnz, ME, lde, m, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, lde, m
    integer, intent(in) :: rowind(1:*), colind(1:*)
    double precision, intent(in) :: spMA(1:*)
    double precision, intent(out) :: ME(1:lde,1:n)
    integer, intent(out) :: info

    double precision, allocatable :: MA(:,:)

    allocate(MA(1:n,1:n))

    call coo_to_dense(n, n, spMA, rowind, colind, nnz, MA, n)
    call mpow_dense(trans, n, MA, n, ME, lde, m, info)

    deallocate(MA)
  end subroutine mpow_coo

end module mpow

