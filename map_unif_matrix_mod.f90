!
! unif matrix
!

module map_unif_matrix
  implicit none

contains

  !
  ! Description:
  !
  !   Generate uniformized kernel from D0, D1
  !
  ! Parameters:
  !   n: size of square matrix A
  !

  subroutine map_unif_dense(n, D0, ldd0, D1, ldd1, P0, ldp0, P1, ldp1, qv, ufact)
    integer, intent(in) :: n, ldd0, ldd1, ldp0, ldp1
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: D0(1:ldd0,1:n), D1(1:ldd1,1:n)
    double precision, intent(out) :: P0(1:ldp0,1:n), P1(1:ldp1,1:n)
    integer :: i

    qv = 0.0d0
    do i = 1, n
      if (qv < -D0(i,i)) then
        qv = -D0(i,i)
      end if
    end do
    qv = qv * ufact

    P0(1:n,1:n) = D0(1:n,1:n)
    P1(1:n,1:n) = D1(1:n,1:n)
    do i = 1, n
      call dscal(n, 1.0d0/qv, P0(1,i), 1)
      call dscal(n, 1.0d0/qv, P1(1,i), 1)
    end do
    do i = 1, n
      P0(i,i) = P0(i,i) + 1.0d0
    end do
  end subroutine map_unif_dense

  subroutine map_unif_csr(n, D0, rowptr0, colind0, nnz0, &
    D1, rowptr1, colind1, nnz1, P0, P1, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: D0(base:base+nnz0-1), D1(base:base+nnz1-1)
    double precision, intent(out) :: P0(base:base+nnz0-1), P1(base:base+nnz1-1)
    integer, intent(in) :: rowptr0(base:*), colind0(base:*)
    integer, intent(in) :: rowptr1(base:*), colind1(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do i = base, base+n-1
      do z = rowptr0(i), rowptr0(i+1)-1
        j = colind0(z)
        if (i == j) then
          if (qv < -D0(z)) then
            qv = -D0(z)
          end if
          diag(i) = z
          exit
        end if
      end do
    end do
    qv = qv * ufact

    P0 = D0
    P1 = D1
    call dscal(nnz0, 1.0d0/qv, P0, 1)
    call dscal(nnz1, 1.0d0/qv, P1, 1)
    do i = base, base+n-1
      P0(diag(i)) = P0(diag(i)) + 1.0d0
    end do
  end subroutine map_unif_csr

  subroutine map_unif_csc(n, D0, colptr0, rowind0, nnz0, &
    D1, colptr1, rowind1, nnz1, P0, P1, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: D0(base:base+nnz0-1), D1(base:base+nnz1-1)
    double precision, intent(out) :: P0(base:base+nnz0-1), P1(base:base+nnz1-1)
    integer, intent(in) :: colptr0(base:*), rowind0(base:*)
    integer, intent(in) :: colptr1(base:*), rowind1(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do j = base, base+n-1
      do z = colptr0(j), colptr0(j+1)-1
        i = rowind0(z)
        if (i == j) then
          if (qv < -D0(z)) then
            qv = -D0(z)
          end if
          diag(i) = z
          exit
        end if
      end do
    end do
    qv = qv * ufact

    P0 = D0
    P1 = D1
    call dscal(nnz0, 1.0d0/qv, P0, 1)
    call dscal(nnz1, 1.0d0/qv, P1, 1)
    do i = base, base+n-1
      P0(diag(i)) = P0(diag(i)) + 1.0d0
    end do
  end subroutine map_unif_csc

  subroutine map_unif_coo(n, D0, rowind0, colind0, nnz0, &
    D1, rowind1, colind1, nnz1, P0, P1, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz0, nnz1
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: D0(base:base+nnz0-1), D1(base:base+nnz1-1)
    double precision, intent(out) :: P0(base:base+nnz0-1), P1(base:base+nnz1-1)
    integer, intent(in) :: rowind0(base:*), colind0(base:*)
    integer, intent(in) :: rowind1(base:*), colind1(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do z = base, base+nnz0-1
      i = rowind0(z)
      j = colind0(z)
      if (i == j) then
        if (qv < -D0(z)) then
          qv = -D0(z)
        end if
        diag(i) = z
      end if
    end do
    qv = qv * ufact

    P0 = D0
    P1 = D1
    call dscal(nnz0, 1.0d0/qv, P0, 1)
    call dscal(nnz1, 1.0d0/qv, P1, 1)
    do i = base, base+n-1
      P0(diag(i)) = P0(diag(i)) + 1.0d0
    end do
  end subroutine map_unif_coo

end module map_unif_matrix

