!
! sparse matrix module
!

module sparse
  implicit none
  logical, public, parameter :: Cindex = .true.
  integer, public, parameter :: sparse_base_index = 0 ! Cindex == .true. -> 0, otherwise -> 1

contains

  ! Description:
  !   return the number of non-zero elements (nnz) of matrix A
  ! Parameters:
  !   m: rows of matrix A
  !   n: columns of matrix A
  !   A: 2D-array of A

  function sparse_nnz(m, n, A, lda) result (nnz)
    integer, intent(in) :: m, n, lda
    double precision, intent(in) :: A(1:lda,1:n)
    integer :: i, j, nnz
    nnz = 0
    do j = 1, n
      do i = 1, m
        if (A(i,j) /= 0.0d0) then
          nnz = nnz + 1
        end if
      end do
    end do
  end function sparse_nnz

  ! Description: Transrate from a dense matrix A to CSR format
  ! Parameters:
  !    m: rows of A
  !    n: columns of A
  !  nnz: nnz of A
  !    A: dense matrix
  !  rowptr: rowptr of csr
  !  colind: colind of csr
  !  spA: value

  subroutine dense_to_csr(m, n, A, lda, spA, rowptr, colind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: A(base:base+lda-1,base:base+n-1)
    double precision, intent(out) :: spA(base:base+nnz-1)
    integer, intent(out) :: rowptr(base:base+m), colind(base:base+nnz-1)

    integer :: i, j, z

    z = base
    rowptr(base) = z
    do i = base, base+m-1
       do j = base, base+n-1
          if (A(i,j) /= 0.0d0) then
             spA(z) = A(i,j)
             colind(z) = j
             z = z + 1
          end if
       end do
       rowptr(i+1) = z
    end do
  end subroutine dense_to_csr

  ! Description: Transrate from a dense matrix A to CSC format
  ! Parameters:
  !    m: rows of A
  !    n: columns of A
  !  nnz: nnz of A
  !    A: dense matrix
  !  colptr: colptr of csc
  !  rowind: rowind of csc
  !  spA: value

  subroutine dense_to_csc(m, n, A, lda, spA, colptr, rowind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: A(base:base+lda-1,base:base+n-1)
    double precision, intent(out) :: spA(base:base+nnz-1)
    integer, intent(out) :: colptr(base:base+n), rowind(base:base+nnz-1)

    integer :: i, j, z

    z = base
    colptr(base) = z
    do j = base, base+n-1
       do i = base, base+m-1
          if (A(i,j) /= 0.0d0) then
             spA(z) = A(i,j)
             rowind(z) = i
             z = z + 1
          end if
       end do
       colptr(j+1) = z
    end do
  end subroutine dense_to_csc

  ! Description: Transrate from a dense matrix A to COO format
  ! Parameters:
  !    m: rows of A
  !    n: columns of A
  !  nnz: nnz of A
  !    A: dense matrix
  !  rowind: colptr of coo
  !  colind: rowind of coo
  !  spA: value

  subroutine dense_to_coo(m, n, A, lda, spA, rowind, colind, nnz)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: A(base:base+lda-1,base:base+n-1)
    double precision, intent(out) :: spA(base:base+nnz-1)
    integer, intent(out) :: rowind(base:base+nnz-1), colind(base:base+nnz-1)

    integer :: i, j, z

    z = base
    do j = base, base+n-1
       do i = base, base+m-1
          if (A(i,j) /= 0.0d0) then
             spA(z) = A(i,j)
             rowind(z) = i
             colind(z) = j
             z = z + 1
          end if
       end do
    end do
  end subroutine dense_to_coo

  ! Description: Transrate from a spase matrix to a
  !              dense matrix
  ! Parameters:
  !    m: rows of spA
  !    n: columns of spA
  !    nnz: nnz of spA
  !    spA: value of spA
  !    rowptr: rowptr of spA
  !    colind: colind of spA
  !    A: dense matrix
  !

  subroutine csr_to_dense(m, n, spA, rowptr, colind, nnz, A, lda)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: spA(base:base+nnz-1)
    integer, intent(in) :: rowptr(base:base+m), colind(base:base+nnz-1)
    double precision, intent(out) :: A(base:base+lda-1,base:base+n-1)

    integer :: i, z

    A(base:base+m-1,base:base+n-1) = 0.0d0
    do i = base, base+m-1
       do z = rowptr(i), rowptr(i+1)-1
          A(i,colind(z)) = spA(z)
       end do
    end do
  end subroutine csr_to_dense

  ! Description: Transrate from a spase matrix to a
  !              dense matrix
  ! Parameters:
  !    m: rows of spA
  !    n: columns of spA
  !    nnz: nnz of spA
  !    spA: value of spA
  !    colptr: rowptr of spA
  !    rowind: colind of spA
  !    A: dense matrix
  !

  subroutine csc_to_dense(m, n, spA, colptr, rowind, nnz, A, lda)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: spA(base:base+nnz-1)
    integer, intent(in) :: colptr(base:base+n), rowind(base:base+nnz-1)
    double precision, intent(out) :: A(base:base+lda-1,base:base+n-1)

    integer :: j, z

    A(base:base+m-1,base:base+n-1) = 0.0d0
    do j = base, base+n-1
       do z = colptr(j), colptr(j+1)-1
          A(rowind(z),j) = spA(z)
       end do
    end do
  end subroutine csc_to_dense

  ! Description: Transrate from a spase matrix to a
  !              dense matrix
  ! Parameters:
  !    m: rows of spA
  !    n: columns of spA
  !    nnz: nnz of spA
  !    spA: value of spA
  !    rowind: rowptr of spA
  !    colind: colind of spA
  !    A: dense matrix
  !

  subroutine coo_to_dense(m, n, spA, rowind, colind, nnz, A, lda)
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: m, n, nnz, lda
    double precision, intent(in) :: spA(base:base+nnz-1)
    integer, intent(in) :: rowind(base:base+nnz-1), colind(base:base+nnz-1)
    double precision, intent(out) :: A(base:base+lda-1,base:base+n-1)

    integer :: z

    A(base:base+m-1,base:base+n-1) = 0.0d0
    do z = base, base+nnz-1
      A(rowind(z),colind(z)) = spA(z)
    end do
  end subroutine coo_to_dense

end module sparse

