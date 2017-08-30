!
! dense kron matrix module
!

module kron
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

  subroutine dkronmv(trans, m, n, n1, n2, alpha, &
    A, lda, x, beta, y)

    character, intent(in) :: trans
    integer, intent(in) :: m, n, n1, n2, lda
    double precision, intent(in) :: alpha, beta, A(1:lda,1:n), x(1:*)
    double precision, intent(out) :: y(1:*)
    integer :: i, j, k

    select case (trans)
      case ('N', 'n')
        call dscal(m*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dgemm('N', 'T', n2, m, n, alpha, x(i), n2, A, lda, 1.0d0, y(j), n2)
          i = i + n*n2
          j = j + m*n2
        end do
      case ('T', 't')
        call dscal(n*n1*n2, beta, y, 1)
        i = 1
        j = 1
        do k = 1, n1
          call dgemm('N', 'N', n2, n, m, alpha, x(i), n2, A, lda, 1.0d0, y(j), n2)
          i = i + m*n2
          j = j + n*n2
        end do
      case default
        write (0,*), 'Unknown character for trans.'
    end select
  end subroutine dkronmv

  !
  ! Description:
  !   This routine computes the following aggregation on kronecker product:
  ! 
  !     B = B + alpha x^T y
  !       and element-wise aggregation of the m-by-n matrix A; I_{n1} * A * I_{n2}
  !
  !   where A is m-by-n matrix
  !         I_{n} is an n-by-n identity matrix, and;
  !         * is kronecker product.
  !
  ! Parameters:
  !   m, n: row and column sizes of square matrix A
  !   n1: size of identity matrix located on left
  !   n2: size of identity matrix located on right
  !

  subroutine dkronr(m, n, n1, n2, alpha, x, y, A, lda)
    integer, intent(in) :: m, n, n1, n2, lda
    double precision, intent(in) :: alpha, x(1:*), y(1:*)
    double precision, intent(out) :: A(1:lda,1:n)
    integer :: i, j, k

    i = 1
    j = 1
    do k = 1, n1
      call dgemm('T', 'N', m, n, n2, alpha, x(i), n2, y(j), n2, 1.0, A, lda)
      i = i + m*n2
      j = j + n*n2
    end do
  end subroutine dkronr

end module kron

