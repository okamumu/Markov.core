!
! sparse kron matrix module
!

module spkron
  use spblas
  implicit none

contains

  !
  ! Description:
  !   This routine computes the following kronecker product:
  ! 
  !     y = alpha (I_{n1} * trans(A) * I_{n2}) x + beta y
  !
  !   where A is an m-by-n square matrix,
  !         I_{n} is an n-by-n identity matrix, and;
  !         * is kronecker product.
  !
  ! Parameters:
  !   m, n: rows and columns of matrix A
  !   n1: size of identity matrix located on left
  !   n2: size of identity matrix located on right
  !

  subroutine dcsrkronmv(trans, m, n, n1, n2, alpha, &
    spA, rowptr, colind, nnz, x, beta, y)

    character, intent(in) :: trans
    integer, intent(in) :: m, n, n1, n2, nnz
    double precision, intent(in) :: alpha, beta, spA(1:*), x(1:*)
    double precision, intent(out) :: y(1:*)
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    integer :: i, j, k

    select case (trans)
      case ('N', 'n')
        call dscal(m*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcsrmm2('N', 'T', n2, m, n, alpha, x(i), n2, spA, rowptr, colind, nnz, 1.0d0, y(j), n2)
          i = i + n*n2
          j = j + m*n2
        end do
      case ('T', 't')
        call dscal(n*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcsrmm2('N', 'N', n2, n, m, alpha, x(i), n2, spA, rowptr, colind, nnz, 1.0d0, y(j), n2)
          i = i + m*n2
          j = j + n*n2
        end do
      case default
        write (0,*), 'Unknown character for trans.'
    end select
  end subroutine dcsrkronmv

  subroutine dcsckronmv(trans, m, n, n1, n2, alpha, &
    spA, colptr, rowind, nnz, x, beta, y)

    character, intent(in) :: trans
    integer, intent(in) :: m, n, n1, n2, nnz
    double precision, intent(in) :: alpha, beta, spA(1:*), x(1:*)
    double precision, intent(out) :: y(1:*)
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    integer :: i, j, k

    select case (trans)
      case ('N', 'n')
        call dscal(m*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcscmm2('N', 'T', n2, m, n, alpha, x(i), n2, spA, colptr, rowind, nnz, 1.0d0, y(j), n2)
          i = i + n*n2
          j = j + m*n2
        end do
      case ('T', 't')
        call dscal(n*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcscmm2('N', 'N', n2, n, m, alpha, x(i), n2, spA, colptr, rowind, nnz, 1.0d0, y(j), n2)
          i = i + m*n2
          j = j + n*n2
        end do
      case default
        write (0,*), 'Unknown character for trans.'
    end select
  end subroutine dcsckronmv

  subroutine dcookronmv(trans, m, n, n1, n2, alpha, &
    spA, rowind, colind, nnz, x, beta, y)

    character, intent(in) :: trans
    integer, intent(in) :: m, n, n1, n2, nnz
    double precision, intent(in) :: alpha, beta, spA(1:*), x(1:*)
    double precision, intent(out) :: y(1:*)
    integer, intent(in) :: rowind(1:*), colind(1:*)
    integer :: i, j, k

    select case (trans)
      case ('N', 'n')
        call dscal(m*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcoomm2('N', 'T', n2, m, n, alpha, x(i), n2, spA, rowind, colind, nnz, 1.0d0, y(j), n2)
          i = i + n*n2
          j = j + m*n2
        end do
      case ('T', 't')
        call dscal(n*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dcoomm2('N', 'N', n2, n, m, alpha, x(i), n2, spA, rowind, colind, nnz, 1.0d0, y(j), n2)
          i = i + m*n2
          j = j + n*n2
        end do
      case default
        write (0,*), 'Unknown character for trans.'
    end select
  end subroutine dcookronmv

end module spkron

