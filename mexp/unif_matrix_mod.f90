!
! unif matrix
!

module unif_matrix
  implicit none

contains

  !
  ! Description:
  !
  !   Generate uniformized kernel from Q
  !
  ! Parameters:
  !   n: size of square matrix A
  !

  subroutine unif_dense(n, Q, ldq, P, ldp, qv, ufact)
    integer, intent(in) :: n, ldq, ldp
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: Q(1:ldq,1:n)
    double precision, intent(out) :: P(1:ldp,1:n)
    integer :: i

    qv = 0.0d0
    do i = 1, n
      if (qv < -Q(i,i)) then
        qv = -Q(i,i)
      end if
    end do
    qv = qv * ufact
    P(1:n,1:n) = Q(1:n,1:n)
    do i = 1, n
      call dscal(n, 1.0d0/qv, P(1,i), 1)
    end do
    do i = 1, n
      P(i,i) = P(i,i) + 1.0d0
    end do
  end subroutine unif_dense

  subroutine unif_csr(n, spQ, rowptr, colind, nnz, spP, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: spQ(base:base+nnz-1)
    double precision, intent(out) :: spP(base:base+nnz-1)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do i = base, base+n-1
      do z = rowptr(i), rowptr(i+1)-1
        j = colind(z)
        if (i == j) then
          if (qv < -spQ(z)) then
            qv = -spQ(z)
          end if
          diag(i) = z
          exit
        end if
      end do
    end do
    qv = qv * ufact

    spP = spQ
    call dscal(nnz, 1.0d0/qv, spP, 1)
    do i = base, base+n-1
      spP(diag(i)) = spP(diag(i)) + 1.0d0
    end do
  end subroutine unif_csr

  subroutine unif_csc(n, spQ, colptr, rowind, nnz, spP, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: spQ(base:base+nnz-1)
    double precision, intent(out) :: spP(base:base+nnz-1)
    integer, intent(in) :: colptr(base:*), rowind(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do j = base, base+n-1
      do z = colptr(j), colptr(j+1)-1
        i = rowind(z)
        if (i == j) then
          if (qv < -spQ(z)) then
            qv = -spQ(z)
          end if
          diag(j) = z
          exit
        end if
      end do
    end do
    qv = qv * ufact

    spP = spQ
    call dscal(nnz, 1.0d0/qv, spP, 1)
    do j = base, base+n-1
      spP(diag(j)) = spP(diag(j)) + 1.0d0
    end do
  end subroutine unif_csc

  subroutine unif_coo(n, spQ, rowind, colind, nnz, spP, qv, ufact)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz
    double precision, intent(in) :: ufact
    double precision, intent(out) :: qv
    double precision, intent(in) :: spQ(base:base+nnz-1)
    double precision, intent(out) :: spP(base:base+nnz-1)
    integer, intent(in) :: rowind(base:*), colind(base:*)
    integer :: i, j, z, diag(base:base+n-1)

    qv = 0.0d0
    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      if (i == j) then
        if (qv < -spQ(z)) then
          qv = -spQ(z)
        end if
        diag(i) = z
      end if
    end do
    qv = qv * ufact

    spP = spQ
    call dscal(nnz, 1.0d0/qv, spP, 1)
    do i = base, base+n-1
      spP(diag(i)) = spP(diag(i)) + 1.0d0
    end do
  end subroutine unif_coo

end module unif_matrix

