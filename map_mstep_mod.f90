!
! gph mstep
!

module map_mstep
  implicit none

contains

  subroutine map_mstep_dense(n, eb, ez, en0, en1, alpha, D0, D1)
    integer, intent(in) :: n
    double precision, intent(in) :: eb(1:n), ez(1:n)
    double precision, intent(in) :: en0(1:n,1:n), en1(1:n,1:n)
    double precision, intent(out) :: alpha(1:n), D0(1:n,1:n), D1(1:n,1:n)
    integer :: i, j
    double precision :: tmp(1:n)

    alpha(1:n) = eb(1:n)
    tmp(1:n) = 0.0d0
    do j = 1, n
      do i = 1, n
        if (i /= j) then
          D0(i,j) = en0(i,j) / ez(i)
          tmp(i) = tmp(i) + D0(i,j)
        end if
        D1(i,j) = en1(i,j) / ez(i)
        tmp(i) = tmp(i) + D1(i,j)
      end do
    end do

    do i = 1, n
      D0(i,i) = -tmp(i)
    end do
  end subroutine map_mstep_dense

  subroutine map_mstep_csr(n, eb, ez, en0, rowptr0, colind0, nnz0, &
    en1, rowptr1, colind1, nnz1, alpha, D0, D1)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en0(base:base+nnz0-1)
    double precision, intent(in) :: en1(base:base+nnz1-1)
    integer, intent(in) :: rowptr0(base:*), colind0(base:*)
    integer, intent(in) :: rowptr1(base:*), colind1(base:*)
    double precision, intent(out) :: alpha(1:n)
    double precision, intent(out) :: D0(base:base+nnz0-1)
    double precision, intent(out) :: D1(base:base+nnz1-1)

    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:)
    tmp(:) = 0.0d0

    do i = base, base+n-1
      do z = rowptr0(i), rowptr0(i+1)-1
        j = colind0(z)
        if (i /= j) then
          D0(z) = en0(z) / ez(i)
          tmp(i) = tmp(i) + D0(z)
        else
          diag(i) = z
        end if
      end do
      do z = rowptr1(i), rowptr1(i+1)-1
        j = colind1(z)
        D1(z) = en1(z) / ez(i)
        tmp(i) = tmp(i) + D1(z)
      end do
    end do

    do i = base, base+n-1
      D0(diag(i)) = -tmp(i)
    end do
  end subroutine map_mstep_csr

  subroutine map_mstep_csc(n, eb, ez, en0, colptr0, rowind0, nnz0, &
    en1, colptr1, rowind1, nnz1, alpha, D0, D1)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en0(base:base+nnz0-1)
    double precision, intent(in) :: en1(base:base+nnz1-1)
    integer, intent(in) :: colptr0(base:*), rowind0(base:*)
    integer, intent(in) :: colptr1(base:*), rowind1(base:*)
    double precision, intent(out) :: alpha(1:n)
    double precision, intent(out) :: D0(base:base+nnz0-1)
    double precision, intent(out) :: D1(base:base+nnz1-1)

    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:)
    tmp(:) = 0.0d0

    do j = base, base+n-1
      do z = colptr0(j), colptr0(j+1)-1
        i = rowind0(z)
        if (i /= j) then
          D0(z) = en0(z) / ez(i)
          tmp(i) = tmp(i) + D0(z)
        else
          diag(i) = z
        end if
      end do
      do z = colptr1(j), colptr1(j+1)-1
        i = rowind1(z)
        D1(z) = en1(z) / ez(i)
        tmp(i) = tmp(i) + D1(z)
      end do
    end do

    do i = base, base+n-1
      D0(diag(i)) = -tmp(i)
    end do
  end subroutine map_mstep_csc

  subroutine map_mstep_coo(n, eb, ez, en0, rowind0, colind0, nnz0, &
    en1, rowind1, colind1, nnz1, alpha, D0, D1)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en0(base:base+nnz0-1)
    double precision, intent(in) :: en1(base:base+nnz1-1)
    integer, intent(in) :: rowind0(base:*), colind0(base:*)
    integer, intent(in) :: rowind1(base:*), colind1(base:*)
    double precision, intent(out) :: alpha(1:n)
    double precision, intent(out) :: D0(base:base+nnz0-1)
    double precision, intent(out) :: D1(base:base+nnz1-1)

    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:)
    tmp(:) = 0.0d0

    do z = base, base+nnz0-1
      i = rowind0(z)
      j = colind0(z)
      if (i /= j) then
        D0(z) = en0(z) / ez(i)
        tmp(i) = tmp(i) + D0(z)
      else
        diag(i) = z
      end if
    end do
    do z = base, base+nnz1-1
      i = rowind1(z)
      j = colind1(z)
      D1(z) = en1(z) / ez(i)
      tmp(i) = tmp(i) + D1(z)
    end do

    do i = base, base+n-1
      D0(diag(i)) = -tmp(i)
    end do
  end subroutine map_mstep_coo

end module map_mstep

