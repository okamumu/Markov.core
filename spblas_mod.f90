!
! sparse BLAS
!

module spblas
  use sparse
  implicit none

contains

  ! Description:
  !     vector-matrix operation for sparsematrix;
  !           y = alpha * trans(spA) * x + beta * y
  ! Parameters:
  !     m: rows of spA
  !     n: columns of spA
  !     nnz: nnz of spA
  !     spA: value of spA
  !     rowptr: rowptr of spA
  !     colptr: colptr of spA
  !     rowind: rowind of spA
  !     colind: colind of spA
  !     alpha, beta: constant values
  !     x, y: vectors
  !

  subroutine dcsrmv(trans, m, n, alpha, &
    spA, rowptr, colind, nnz, x, incx, beta, y, incy)
    integer, parameter :: base = sparse_base_index
    character, intent(in) :: trans
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, beta, spA(base:*), x(1:incx,base:*)
    double precision, intent(out) :: y(1:incy,base:*)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    integer :: i, j, z

    select case (trans)
      case ('N', 'n')
        call dscal(m, beta, y(1,base), incy)
        do i = base, base+m-1
           do z = rowptr(i), rowptr(i+1)-1
              j = colind(z)
              y(1,i) = y(1,i) + alpha * spA(z) * x(1,j)
           end do
        end do
      case ('T', 't')
        call dscal(n, beta, y(1,base), incy)
        do i = base, base+m-1
           do z = rowptr(i), rowptr(i+1)-1
              j = colind(z)
              y(1,j) = y(1,j) + alpha * spA(z) * x(1,i)
           end do
        end do
    end select
  end subroutine dcsrmv

  subroutine dcscmv(trans, m, n, alpha, &
    spA, colptr, rowind, nnz, x, incx, beta, y, incy)
    integer, parameter :: base = sparse_base_index
    character, intent(in) :: trans
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, beta, spA(base:*), x(1:incx,base:*)
    double precision, intent(out) :: y(1:incy,base:*)
    integer, intent(in) :: colptr(base:*), rowind(base:*)
    integer :: i, j, z

    select case (trans)
      case ('N', 'n')
        call dscal(m, beta, y(1,base), incy)
        do j = base, base+n-1
           do z = colptr(j), colptr(j+1)-1
              i = rowind(z)
              y(1,i) = y(1,i) + alpha * spA(z) * x(1,j)
           end do
        end do
      case ('T', 't')
        call dscal(n, beta, y(1,base), incy)
        do j = base, base+n-1
           do z = colptr(j), colptr(j+1)-1
              i = rowind(z)
              y(1,j) = y(1,j) + alpha * spA(z) * x(1,i)
           end do
        end do
    end select
  end subroutine dcscmv

  subroutine dcoomv(trans, m, n, alpha, &
    spA, rowind, colind, nnz, x, incx, beta, y, incy)
    integer, parameter :: base = sparse_base_index
    character, intent(in) :: trans
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, beta, spA(base:*), x(1:incx,base:*)
    double precision, intent(out) :: y(1:incy,base:*)
    integer, intent(in) :: rowind(base:*), colind(base:*)
    integer :: i, j, z

    select case (trans)
      case ('N', 'n')
        call dscal(m, beta, y(1,base), incy)
        do z = base, base+nnz-1
          i = rowind(z)
          j = colind(z)
          y(1,i) = y(1,i) + alpha * spA(z) * x(1,j)
        end do
      case ('T', 't')
        call dscal(n, beta, y(1,base), incy)
        do z = base, base+nnz-1
          i = rowind(z)
          j = colind(z)
          y(1,j) = y(1,j) + alpha * spA(z) * x(1,i)
        end do
    end select
  end subroutine dcoomv

  ! Description:
  !     vector-vector operation for sparsematrix;
  !           spA = alpha * x * y + spA
  ! Parameters:
  !     m: rows of spA
  !     n: columns of spA
  !     nnz: nnz of spA
  !     spA: value of spA
  !     rowptr: rowptr of spA
  !     colptr: colptr of spA
  !     rowind: rowind of spA
  !     colind: colind of spA
  !     alpha, beta: constant values
  !     x, y: vectors
  !

  subroutine dcsrr(m, n, alpha, x, incx, y, incy, &
    spA, rowptr, colind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, x(1:incx,base:*), y(1:incy,base:*)
    double precision, intent(out) :: spA(base:*)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    integer :: i, j, z
    do i = base, base+m-1
       do z = rowptr(i), rowptr(i+1)-1
          j = colind(z)
          spA(z) = spA(z) + alpha * x(1,i) * y(1,j)
       end do
    end do
  end subroutine dcsrr

  subroutine dcscr(m, n, alpha, x, incx, y, incy, &
    spA, colptr, rowind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, x(1:incx,base:*), y(1:incy,base:*)
    double precision, intent(out) :: spA(base:*)
    integer, intent(in) :: colptr(base:*), rowind(base:*)
    integer :: i, j, z
    do j = base, base+n-1
       do z = colptr(j), colptr(j+1)-1
          i = rowind(z)
          spA(z) = spA(z) + alpha * x(1,i) * y(1,j)
       end do
    end do
  end subroutine dcscr

  subroutine dcoor(m, n, alpha, x, incx, y, incy, &
    spA, rowind, colind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, incx, incy
    double precision, intent(in) :: alpha, x(1:incx,base:*), y(1:incy,base:*)
    double precision, intent(out) :: spA(base:*)
    integer, intent(in) :: rowind(base:*), colind(base:*)
    integer :: i, j, z
    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      spA(z) = spA(z) + alpha * x(1,i) * y(1,j)
    end do
  end subroutine dcoor

  ! Description:
  !     matrix-matrix operation for sparsematrix;
  !           C = alpha * transA(spA) * transB(B) * x + beta * C
  ! Parameters:
  !     m: rows of spA
  !     n: columns of spA
  !     nnz: nnz of spA
  !     spA: value of spA
  !     rowptr: rowptr of spA
  !     colptr: colptr of spA
  !     rowind: rowind of spA
  !     colind: colind of spA
  !     alpha, beta: constant values
  !     x, y: vectors
  !

  subroutine dcsrmm(transA, transB, m, n, k, alpha, &
    spA, rowptr, colind, nnz, B, ldb, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, ldb, ldc
    double precision, intent(in) :: alpha, beta, spA(1:*), B(1:ldb,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    integer :: i

    select case (transB)
      case ('N', 'n')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcsrmv('N', m, k, alpha, spA, rowptr, colind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcsrmv('T', k, m, alpha, spA, rowptr, colind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
        end select
      case ('T', 't')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcsrmv('N', m, k, alpha, spA, rowptr, colind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcsrmv('T', k, m, alpha, spA, rowptr, colind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
        end select
    end select
  end subroutine dcsrmm

  subroutine dcscmm(transA, transB, m, n, k, alpha, &
    spA, colptr, rowind, nnz, B, ldb, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, ldb, ldc
    double precision, intent(in) :: alpha, beta, spA(1:*), B(1:ldb,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    integer :: i

    select case (transB)
      case ('N', 'n')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcscmv('N', m, k, alpha, spA, colptr, rowind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcscmv('T', k, m, alpha, spA, colptr, rowind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
        end select
      case ('T', 't')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcscmv('N', m, k, alpha, spA, colptr, rowind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcscmv('T', k, m, alpha, spA, colptr, rowind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
        end select
    end select
  end subroutine dcscmm

  subroutine dcoomm(transA, transB, m, n, k, alpha, &
    spA, rowind, colind, nnz, B, ldb, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, ldb, ldc
    double precision, intent(in) :: alpha, beta, spA(1:*), B(1:ldb,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: rowind(1:*), colind(1:*)
    integer :: i

    select case (transB)
      case ('N', 'n')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcoomv('N', m, k, alpha, spA, rowind, colind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcoomv('T', k, m, alpha, spA, rowind, colind, nnz, B(1,i), 1, beta, C(1,i), 1)
            end do
        end select
      case ('T', 't')
        select case (transA)
          case ('N', 'n')
            do i = 1, n
              call dcoomv('N', m, k, alpha, spA, rowind, colind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
          case ('T', 't')
            do i = 1, n
              call dcoomv('T', k, m, alpha, spA, rowind, colind, nnz, B(i,1), ldb, beta, C(1,i), 1)
            end do
        end select
    end select
  end subroutine dcoomm

  ! Description:
  !     matrix-matrix operation for sparsematrix;
  !           C = alpha * transA(A) * transB(spB) * x + beta * C
  ! Parameters:
  !     m: rows of C
  !     n: columns of C
  !     nnz: nnz of spA
  !     spA: value of spA
  !     rowptr: rowptr of spA
  !     colptr: colptr of spA
  !     rowind: rowind of spA
  !     colind: colind of spA
  !     alpha, beta: constant values
  !     x, y: vectors
  !

  subroutine dcsrmm2(transA, transB, m, n, k, alpha, &
    A, lda, spB, rowptr, colind, nnz, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, lda, ldc
    double precision, intent(in) :: alpha, beta, spB(1:*), A(1:lda,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: rowptr(1:*), colind(1:*)
    integer :: i

    select case (transA)
      case ('N', 'n')
        select case (transB)
          case ('N', 'n')
            do i = 1, m
              call dcsrmv('T', k, n, alpha, spB, rowptr, colind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcsrmv('N', n, k, alpha, spB, rowptr, colind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
        end select
      case ('T', 't')
        select case (transB)
          case ('N', 'n')
            do i = 1, n
              call dcsrmv('T', k, n, alpha, spB, rowptr, colind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcsrmv('N', n, k, alpha, spB, rowptr, colind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
        end select
    end select
  end subroutine dcsrmm2

  subroutine dcscmm2(transA, transB, m, n, k, alpha, &
    A, lda, spB, colptr, rowind, nnz, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, lda, ldc
    double precision, intent(in) :: alpha, beta, spB(1:*), A(1:lda,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: colptr(1:*), rowind(1:*)
    integer :: i

    select case (transA)
      case ('N', 'n')
        select case (transB)
          case ('N', 'n')
            do i = 1, m
              call dcscmv('T', k, n, alpha, spB, colptr, rowind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcscmv('N', n, k, alpha, spB, colptr, rowind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
        end select
      case ('T', 't')
        select case (transB)
          case ('N', 'n')
            do i = 1, n
              call dcscmv('T', k, n, alpha, spB, colptr, rowind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcscmv('N', n, k, alpha, spB, colptr, rowind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
        end select
    end select
  end subroutine dcscmm2

  subroutine dcoomm2(transA, transB, m, n, k, alpha, &
    A, lda, spB, rowind, colind, nnz, beta, C, ldc)

    character, intent(in) :: transA, transB
    integer, intent(in) :: m, n, k, nnz, lda, ldc
    double precision, intent(in) :: alpha, beta, spB(1:*), A(1:lda,1:*)
    double precision, intent(out) :: C(1:ldc,1:*)
    integer, intent(in) :: rowind(1:*), colind(1:*)
    integer :: i

    select case (transA)
      case ('N', 'n')
        select case (transB)
          case ('N', 'n')
            do i = 1, m
              call dcoomv('T', k, n, alpha, spB, rowind, colind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcoomv('N', n, k, alpha, spB, rowind, colind, nnz, A(i,1), lda, beta, C(i,1), ldc)
            end do
        end select
      case ('T', 't')
        select case (transB)
          case ('N', 'n')
            do i = 1, n
              call dcoomv('T', k, n, alpha, spB, rowind, colind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
          case ('T', 't')
            do i = 1, n
              call dcoomv('N', n, k, alpha, spB, rowind, colind, nnz, A(1,i), 1, beta, C(i,1), ldc)
            end do
        end select
    end select
  end subroutine dcoomm2
end module spblas

