
module arnoldi
  implicit none

contains

!     Description:

!     Arnoldi process. The matrix A is approximated by h on the subspace

!       h = v^T A v on Krylov subspace  {x, A*x, A^2*x, ... A^(m-1)*x}

!     Parameters:

!     A (in): a n-by-n square sparse matrix
!     x (in): a vector with n elements
!     h (out): m-by-m square Hessenburg matrix
!     v (out): n-by-m matrix of orthogonal vectors on Krylov subspace {x, A*x, A^2*x, ... A^(m-1)*x}

!     tol (in): tolerance error for happy breakdown
!     ite (in): the number of times of iteration for Arnoldi process to reduce roundoff errors
!               ite=1 is enough in most cases
!     rnorm (out): 2-norm of vector x (and is used for checking for happy breakdown)
!     info (out): the size of effective dimensions of h when a happy breakdown occurs

!     work: working area. The required size is n.

  subroutine arnoldi_dense(trans, n, A, lda, x, incx, m, H, ldh, V, ldv, &
    beta, rnorm, tol, ite, info)
    character, intent(in) :: trans
    integer, intent(in) :: n, lda, incx, m, ldh, ldv, ite
    double precision, intent(in) :: A(1:lda,1:n), x(1:incx,1:n), tol
    double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
    double precision, intent(out) :: beta, rnorm
    integer, intent(out) :: info

    integer :: j, i, l
    double precision :: dnrm2, ddot
    double precision :: r, tmp(1:n)

    rnorm = sqrt(sum(A(1:n,1:n)**2))
    beta = dnrm2(n, x, incx)
    H(1:m,1:m) = 0.0d0
    V(1:n,1:m) = 0.0d0

    call daxpy(n, 1.0d0/beta, x, incx, V(1,1), 1)
    do j = 1, m
      call dgemv(trans, n, n, 1.0d0, A, lda, V(1,j), 1, 0.0d0, tmp, 1)
      do i = 1, j
        do l = 1, ite
          r = ddot(n, V(1,i), 1, tmp, 1)
          H(i,j) = H(i,j) + r
          call daxpy(n, -r, V(1,i), 1, tmp, 1)
        end do
      end do
      if (j /= m) then
        H(j+1,j) = dnrm2(n, tmp, 1)
        if (H(j+1,j) < rnorm * tol) then
          info = j
          return
        end if
        call daxpy(n, 1.0d0/H(j+1,j), tmp, 1, V(1,j+1), 1)
      end if
    end do
    info = m
  end subroutine arnoldi_dense

  subroutine arnoldi_csr(trans, n, A, rowptr, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
    double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
    double precision, intent(out) :: beta, rnorm
    integer, intent(out) :: info

    integer :: j, i, l
    double precision :: dnrm2, ddot
    double precision :: r, tmp(1:n)

    rnorm = dnrm2(nnz, A, 1)
    beta = dnrm2(n, x, incx)
    H(1:m,1:m) = 0.0d0
    V(1:n,1:m) = 0.0d0

    call daxpy(n, 1.0d0/beta, x, incx, V(1,1), 1)
    do j = 1, m
      call dcsrmv(trans, n, n, 1.0d0, A, rowptr, colind, nnz, V(1,j), 1, 0.0d0, tmp, 1)
      do i = 1, j
        do l = 1, ite
          r = ddot(n, V(1,i), 1, tmp, 1)
          H(i,j) = H(i,j) + r
          call daxpy(n, -r, V(1,i), 1, tmp, 1)
        end do
      end do
      if (j /= m) then
        H(j+1,j) = dnrm2(n, tmp, 1)
        if (H(j+1,j) < rnorm * tol) then
          info = j
          return
        end if
        call daxpy(n, 1.0d0/H(j+1,j), tmp, 1, V(1,j+1), 1)
      end if
    end do
    info = m
  end subroutine arnoldi_csr

  subroutine arnoldi_csc(trans, n, A, colptr, rowind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
    double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
    double precision, intent(out) :: beta, rnorm
    integer, intent(out) :: info

    integer :: j, i, l
    double precision :: dnrm2, ddot
    double precision :: r, tmp(1:n)

    rnorm = dnrm2(nnz, A, 1)
    beta = dnrm2(n, x, incx)
    H(1:m,1:m) = 0.0d0
    V(1:n,1:m) = 0.0d0

    call daxpy(n, 1.0d0/beta, x, incx, V(1,1), 1)
    do j = 1, m
      call dcscmv(trans, n, n, 1.0d0, A, colptr, rowind, nnz, V(1,j), 1, 0.0d0, tmp, 1)
      do i = 1, j
        do l = 1, ite
          r = ddot(n, V(1,i), 1, tmp, 1)
          H(i,j) = H(i,j) + r
          call daxpy(n, -r, V(1,i), 1, tmp, 1)
        end do
      end do
      if (j /= m) then
        H(j+1,j) = dnrm2(n, tmp, 1)
        if (H(j+1,j) < rnorm * tol) then
          info = j
          return
        end if
        call daxpy(n, 1.0d0/H(j+1,j), tmp, 1, V(1,j+1), 1)
      end if
    end do
    info = m
  end subroutine arnoldi_csc

  subroutine arnoldi_coo(trans, n, A, rowind, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
    use spblas
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
    double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
    integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
    double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
    double precision, intent(out) :: beta, rnorm
    integer, intent(out) :: info

    integer :: j, i, l
    double precision :: dnrm2, ddot
    double precision :: r, tmp(1:n)

    rnorm = dnrm2(nnz, A, 1)
    beta = dnrm2(n, x, incx)
    H(1:m,1:m) = 0.0d0
    V(1:n,1:m) = 0.0d0

    call daxpy(n, 1.0d0/beta, x, incx, V(1,1), 1)
    do j = 1, m
      call dcoomv(trans, n, n, 1.0d0, A, rowind, colind, nnz, V(1,j), 1, 0.0d0, tmp, 1)
      do i = 1, j
        do l = 1, ite
          r = ddot(n, V(1,i), 1, tmp, 1)
          H(i,j) = H(i,j) + r
          call daxpy(n, -r, V(1,i), 1, tmp, 1)
        end do
      end do
      if (j /= m) then
        H(j+1,j) = dnrm2(n, tmp, 1)
        if (H(j+1,j) < rnorm * tol) then
          info = j
          return
        end if
        call daxpy(n, 1.0d0/H(j+1,j), tmp, 1, V(1,j+1), 1)
      end if
    end do
    info = m
  end subroutine arnoldi_coo

end module arnoldi

