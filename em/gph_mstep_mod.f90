!
! gph mstep
!

module gph_mstep
  implicit none

contains

  subroutine gph_mstep_dense(n, etotal, eb, ey, ez, en, &
    alpha, tau, T)
    integer, intent(in) :: n
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(in) :: en(1:n,1:n)
    double precision, intent(out) :: alpha(1:n), tau(1:n), T(1:n,1:n)
    integer :: i, j
    double precision :: tmp(1:n)

    alpha(1:n) = eb(1:n) / etotal
    tau(1:n) = ey(1:n) / ez(1:n)

    tmp(1:n) = tau(1:n)
    do j = 1, n
      do i = 1, n
        if (i /= j) then
          T(i,j) = en(i,j) / ez(i)
          tmp(i) = tmp(i) + T(i,j)
        end if
      end do
    end do

    do i = 1, n
      T(i,i) = -tmp(i)
    end do
  end subroutine gph_mstep_dense

  subroutine gph_mstep_csr(n, etotal, eb, ey, ez, en, rowptr, colind, nnz, &
    alpha, tau, spT)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en(base:base+nnz-1)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(base:base+nnz-1)
    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:) / etotal
    tau(:) = ey(:) / ez(:)

    tmp(:) = tau(:)
    do i = base, base+n-1
      do z = rowptr(i), rowptr(i+1)-1
        j = colind(z)
        if (i /= j) then
          spT(z) = en(z) / ez(i)
          tmp(i) = tmp(i) + spT(z)
        else
          diag(i) = z
        end if
      end do
    end do

    do i = base, base+n-1
      spT(diag(i)) = -tmp(i)
    end do
  end subroutine gph_mstep_csr

  subroutine gph_mstep_csc(n, etotal, eb, ey, ez, en, colptr, rowind, nnz, &
    alpha, tau, spT)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en(base:base+nnz-1)
    integer, intent(in) :: colptr(base:*), rowind(base:*)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(base:base+nnz-1)
    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:) / etotal
    tau(:) = ey(:) / ez(:)

    tmp(:) = tau(:)
    do j = base, base+n-1
      do z = colptr(j), colptr(j+1)-1
        i = rowind(z)
        if (i /= j) then
          spT(z) = en(z) / ez(i)
          tmp(i) = tmp(i) + spT(z)
        else
          diag(i) = z
        end if
      end do
    end do

    do i = base, base+n-1
      spT(diag(i)) = -tmp(i)
    end do
  end subroutine gph_mstep_csc

  subroutine gph_mstep_coo(n, etotal, eb, ey, ez, en, rowind, colind, nnz, &
    alpha, tau, spT)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(in) :: en(base:base+nnz-1)
    integer, intent(in) :: rowind(base:*), colind(base:*)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(base:base+nnz-1)
    integer :: i, j, z, diag(base:base+n-1)
    double precision :: tmp(base:base+n-1)

    alpha(:) = eb(:) / etotal
    tau(:) = ey(:) / ez(:)

    tmp(:) = tau(:)
    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      if (i /= j) then
        spT(z) = en(z) / ez(i)
        tmp(i) = tmp(i) + spT(z)
      else
        diag(i) = z
      end if
    end do

    do i = base, base+n-1
      spT(diag(i)) = -tmp(i)
    end do
  end subroutine gph_mstep_coo

end module gph_mstep

