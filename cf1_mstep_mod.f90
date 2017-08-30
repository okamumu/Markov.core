!
! cf1 mstep
!

module cf1_mstep
  implicit none
  private cf1_phase_swap, cf1_phase_sort

contains

  subroutine cf1_phase_swap(i, j, n, alpha, incal, lambda, inclam)
    integer, intent(in) :: i, j, n, incal, inclam
    double precision, intent(inout) :: alpha(1:incal,1:n), lambda(1:inclam,1:n)
    double precision :: w, tmp
    w = lambda(1,j) / lambda(1,i)
    alpha(1,i) = alpha(1,i) + (1.0d0 - w) * alpha(1,j)
    alpha(1,j) = w * alpha(1,j)
    tmp = lambda(1,j)
    lambda(1,j) = lambda(1,i)
    lambda(1,i) = tmp
  end subroutine cf1_phase_swap

  subroutine cf1_phase_sort(n, alpha, incal, lambda, inclam)
    integer, intent(in) :: n, incal, inclam
    double precision, intent(inout) :: alpha(1:incal,1:n), lambda(1:inclam,1:n)
    integer :: i, j

    do i = 1, n-1
      if (lambda(1,i) > lambda(1,i+1)) then
        call cf1_phase_swap(i, i+1, n, alpha, incal, lambda, inclam)
        do j = i, 2, -1
          if (lambda(1,j-1) < lambda(1,j)) then
            exit
          end if
          call cf1_phase_swap(j-1, j, n, alpha, incal, lambda, inclam)
        end do
      end if
    end do
  end subroutine cf1_phase_sort

  !!! reshape

  subroutine cf1_reshape_dense(n, alpha, tau, T)
    integer, intent(in) :: n
    double precision, intent(inout) :: alpha(1:n), tau(1:n), T(1:n,1:n)
    integer :: i
    double precision :: lambda(1:n)

    do i = 1, n-1
      lambda(i) = T(i,i+1)
    end do
    lambda(n) = tau(n)

    call cf1_phase_sort(n, alpha, 1, lambda, 1)

    T(1:n,1:n) = 0.0d0
    tau(1:n) = 0.0d0

    do i = 1, n-1
      T(i,i) = -lambda(i)
      T(i,i+1) = lambda(i)
    end do
    T(n,n) = -lambda(n)
    tau(n) = lambda(n)
  end subroutine cf1_reshape_dense

  subroutine cf1_reshape_csr(n, alpha, tau, spT, rowptr, colind, nnz)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(inout) :: alpha(1:n), tau(base:base+n-1)
    double precision, intent(inout) :: spT(base:base+nnz-1)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    integer :: i, j, z
    double precision :: lambda(base:base+n-1)

    do i = base, base+n-1
      do z = rowptr(i), rowptr(i+1)-1
        j = colind(z)
        if (j == i+1) then
          lambda(i) = spT(z)
          exit
        end if
      end do
    end do
    lambda(base+n-1) = tau(base+n-1)

    call cf1_phase_sort(n, alpha, 1, lambda, 1)

    spT(:) = 0.0d0
    tau(:) = 0.0d0

    do i = base, base+n-1
      do z = rowptr(i), rowptr(i+1)-1
        j = colind(z)
        if (i == j) then
          spT(z) = -lambda(i)
        else if (j == i+1) then
          spT(z) = lambda(i)
        end if
      end do
    end do
    tau(base+n-1) = lambda(base+n-1)
  end subroutine cf1_reshape_csr

  subroutine cf1_reshape_csc(n, alpha, tau, spT, colptr, rowind, nnz)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(inout) :: alpha(1:n), tau(base:base+n-1)
    double precision, intent(inout) :: spT(base:base+nnz-1)
    integer, intent(in) :: colptr(base:*), rowind(base:*)
    integer :: i, j, z
    double precision :: lambda(base:base+n-1)

    do j = base, base+n-1
      do z = colptr(j), colptr(j+1)-1
        i = rowind(z)
        if (j == i+1) then
          lambda(i) = spT(z)
          exit
        end if
      end do
    end do
    lambda(base+n-1) = tau(base+n-1)

    call cf1_phase_sort(n, alpha, 1, lambda, 1)

    spT(:) = 0.0d0
    tau(:) = 0.0d0

    do j = base, base+n-1
      do z = colptr(j), colptr(j+1)-1
        i = rowind(z)
        if (i == j) then
          spT(z) = -lambda(i)
        else if (j == i+1) then
          spT(z) = lambda(i)
        end if
      end do
    end do
    tau(base+n-1) = lambda(base+n-1)
  end subroutine cf1_reshape_csc

  subroutine cf1_reshape_coo(n, alpha, tau, spT, rowind, colind, nnz)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(inout) :: alpha(1:n), tau(base:base+n-1)
    double precision, intent(inout) :: spT(base:base+nnz-1)
    integer, intent(in) :: rowind(base:*), colind(base:*)
    integer :: i, j, z
    double precision :: lambda(base:base+n-1)

    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      if (j == i+1) then
        lambda(i) = spT(z)
      end if
    end do
    lambda(base+n-1) = tau(base+n-1)

    call cf1_phase_sort(n, alpha, 1, lambda, 1)

    spT(:) = 0.0d0
    tau(:) = 0.0d0

    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      if (i == j) then
        spT(z) = -lambda(i)
      else if (j == i+1) then
        spT(z) = lambda(i)
      end if
    end do
    tau(base+n-1) = lambda(base+n-1)
  end subroutine cf1_reshape_coo

  !!! mstep

  subroutine cf1_mstep_dense(n, etotal, eb, ey, ez, en, alpha, tau, T)
    use gph_mstep
    integer, intent(in) :: n
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(in) :: en(1:n,1:n)
    double precision, intent(out) :: alpha(1:n), tau(1:n), T(1:n,1:n)

    call gph_mstep_dense(n, etotal, eb, ey, ez, en, alpha, tau, T)
    call cf1_reshape_dense(n, alpha, tau, T)
  end subroutine cf1_mstep_dense

  subroutine cf1_mstep_csr(n, etotal, eb, ey, ez, en, &
    rowptr, colind, nnz, alpha, tau, spT)
    use gph_mstep
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(in) :: en(1:nnz)
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(1:nnz)

    call gph_mstep_csr(n, etotal, eb, ey, ez, en, rowptr, colind, nnz, alpha, tau, spT)
    call cf1_reshape_csr(n, alpha, tau, spT, rowptr, colind, nnz)
  end subroutine cf1_mstep_csr

  subroutine cf1_mstep_csc(n, etotal, eb, ey, ez, en, &
    colptr, rowind, nnz, alpha, tau, spT)
    use gph_mstep
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(in) :: en(1:nnz)
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(1:nnz)

    call gph_mstep_csc(n, etotal, eb, ey, ez, en, colptr, rowind, nnz, alpha, tau, spT)
    call cf1_reshape_csc(n, alpha, tau, spT, colptr, rowind, nnz)
  end subroutine cf1_mstep_csc

  subroutine cf1_mstep_coo(n, etotal, eb, ey, ez, en, &
    rowind, colind, nnz, alpha, tau, spT)
    use gph_mstep
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(in) :: en(1:nnz)
    integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
    double precision, intent(out) :: alpha(1:n), tau(1:n), spT(1:nnz)

    call gph_mstep_coo(n, etotal, eb, ey, ez, en, rowind, colind, nnz, alpha, tau, spT)
    call cf1_reshape_coo(n, alpha, tau, spT, rowind, colind, nnz)
  end subroutine cf1_mstep_coo

end module cf1_mstep

