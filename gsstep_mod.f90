!
! gsstep
!

module gsstep
  implicit none
  public gsstep_fwd_shift_dense
  public gsstep_bwd_shift_dense
  public gsstep_fwd_shift_csr
  public gsstep_bwd_shift_csr
  public gsstep_fwd_shift_csc
  public gsstep_bwd_shift_csc
!!!
  public gsstep_fwd_shift_dense_mat
  public gsstep_bwd_shift_dense_mat
  public gsstep_fwd_shift_csr_mat
  public gsstep_bwd_shift_csr_mat
  public gsstep_fwd_shift_csc_mat
  public gsstep_bwd_shift_csc_mat
!!!
  private gsstep_fwd_shift_csr_notrans
  private gsstep_fwd_shift_csr_trans
  private gsstep_bwd_shift_csr_notrans
  private gsstep_bwd_shift_csr_trans

contains

! Kahan sum

  subroutine kahan_sum(sum, x, c)
    double precision, intent(in) :: x
    double precision, intent(inout) :: sum, c
    double precision :: t, y
    y = x - c
    t = sum + y
    c = (t - sum) - y
    sum = t
  end subroutine kahan_sum

!  SOR step for solving the following linear equation
!
!         alpha * trans(A - sigma I) * x = b
!
!         The step computes
!
!         fwd
!         x := (D/omega + L)^(-1) (b/alpha - (U - D (1-omega)/omega - sigma I) * x)
!
!         bwd
!         x := (D/omega + U)^(-1) (b/alpha - (L - D (1-omega)/omega - sigma I) * x)
!
!         A: square matrix
!         x: vector (in; initial vector for the step, out; updated vector)
!         b: constant vector

subroutine gsstep_fwd_shift_dense(trans, n, alpha, A, lda, sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, A(1:lda,1:n), b(1:incb,1:n)
  double precision, intent(inout) :: x(1:incx,1:n)

  integer :: i, j
  double precision :: tmpx, tmpd, c

  select case (trans)
    case ('N', 'n')
      do i = 1, n
        tmpd = 0.0d0
        c = 0.0d0
        tmpx = b(1,i) / alpha
        do j = 1, n
          if (i == j) then
            tmpd = A(i,j)
!            tmpx = tmpx + sigma * x(1,j)
            call kahan_sum(tmpx, sigma * x(1,j), c)
          else
!            tmpx = tmpx - A(i,j) * x(1,j)
            call kahan_sum(tmpx, - A(i,j) * x(1,j), c)
          end if
        end do
        x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
      end do
    case ('T', 't')
      do j = 1, n
        tmpd = 0.0d0
        c = 0.0d0
        tmpx = b(1,j) / alpha
        do i = 1, n
          if (i == j) then
            tmpd = A(i,j)
!            tmpx = tmpx + sigma * x(1,i)
            call kahan_sum(tmpx, sigma * x(1,j), c)
          else
!            tmpx = tmpx - A(i,j) * x(1,i)
            call kahan_sum(tmpx, - A(i,j) * x(1,i), c)
          end if
        end do
        x(1,j) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,j)
      end do
  end select
end subroutine gsstep_fwd_shift_dense

!!!!

subroutine gsstep_fwd_shift_csr(trans, n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_fwd_shift_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
        sigma, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_fwd_shift_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
        sigma, omega, b, incb, x, incx)
  end select
end subroutine gsstep_fwd_shift_csr

subroutine gsstep_fwd_shift_csc(trans, n, alpha, spA, colptr, rowind, nnz, &
  sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_fwd_shift_csr_trans(n, alpha, spA, colptr, rowind, nnz, &
        sigma, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_fwd_shift_csr_notrans(n, alpha, spA, colptr, rowind, nnz, &
        sigma, omega, b, incb, x, incx)
  end select
end subroutine gsstep_fwd_shift_csc

!!! private

subroutine gsstep_fwd_shift_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  use sparse
  integer, parameter :: base = sparse_base_index
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(base:*), b(1:incb,base:base+n-1)
  integer, intent(in) :: rowptr(base:*), colind(base:*)
  double precision, intent(inout) :: x(1:incx,base:base+n-1)

  integer :: i, j, z
  double precision :: tmpd, tmpx
  double precision :: c

  do i = base, base+n-1
    tmpd = 0.0d0
    c = 0.0d0
    tmpx = b(1,i) / alpha
    do z = rowptr(i), rowptr(i+1)-1
      j = colind(z)
      if (i == j) then
        tmpd = spA(z)
!!        tmpx = tmpx + sigma * x(1,i)
        call kahan_sum(tmpx, sigma * x(1,i), c)
      else
!!        tmpx = tmpx - spA(z) * x(1,j)
        call kahan_sum(tmpx, - spA(z) * x(1,j), c)
      end if
    end do
    x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
  end do
end subroutine gsstep_fwd_shift_csr_notrans

subroutine gsstep_fwd_shift_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  use sparse
  integer, parameter :: base = sparse_base_index
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(base:*), b(1:incb,base:base+n-1)
  integer, intent(in) :: rowptr(base:*), colind(base:*)
  double precision, intent(inout) :: x(1:incx,base:base+n-1)

  integer :: i, j, z
  double precision :: tmp(base:base+n-1)
  double precision :: c(base:base+n-1)
  integer :: diag(base:base+n-1)

  tmp = 0.0d0
  c = 0.0d0
  call daxpy(n, 1.0/alpha, b, incb, tmp, 1)
  do i = base, base+n-1
    do z = rowptr(i), rowptr(i+1)-1
      j = colind(z)
      if (i == j) then
        diag(i) = z
!        tmp(i) = tmp(i) + sigma * x(1,i)
        call kahan_sum(tmp(i), sigma * x(1,i), c(i))
        exit
      else
!        tmp(j) = tmp(j) - spA(z) * x(1,i)
        call kahan_sum(tmp(j), - spA(z) * x(1,i), c(j))
      end if
    end do
  end do

  do i = base, base+n-1
    x(1,i) = omega * tmp(i) / spA(diag(i)) + (1.0d0 - omega) * x(1,i)
    do z = diag(i)+1, rowptr(i+1)-1
      j = colind(z)
!      tmp(j) = tmp(j) - spA(z) * x(1,i)
      call kahan_sum(tmp(j), - spA(z) * x(1,i), c(j))
    end do
  end do
end subroutine gsstep_fwd_shift_csr_trans

!!!! bwd

subroutine gsstep_bwd_shift_dense(trans, n, alpha, A, lda, sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, A(1:lda,1:n), b(1:incb,1:n)
  double precision, intent(inout) :: x(1:incx,1:n)

  integer :: i, j
  double precision :: tmpx, tmpd, c

  select case (trans)
    case ('N', 'n')
      do i = n, 1, -1
        tmpd = 0.0d0
        c = 0.0d0
        tmpx = b(1,i) / alpha
        do j = n, 1, -1
          if (i == j) then
            tmpd = A(i,j)
!            tmpx = tmpx + sigma * x(1,j)
            call kahan_sum(tmpx, sigma * x(1,j), c)
          else
!            tmpx = tmpx - A(i,j) * x(1,j)
            call kahan_sum(tmpx, - A(i,j) * x(1,j), c)
          end if
        end do
        x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
      end do
    case ('T', 't')
      do j = n, 1, -1
        tmpd = 0.0d0
        c = 0.0d0
        tmpx = b(1,j) / alpha
        do i = n, 1, -1
          if (i == j) then
            tmpd = A(i,j)
!            tmpx = tmpx + sigma * x(1,i)
            call kahan_sum(tmpx, sigma * x(1,i), c)
          else
!            tmpx = tmpx - A(i,j) * x(1,i)
            call kahan_sum(tmpx, - A(i,j) * x(1,i), c)
          end if
        end do
        x(1,j) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,j)
      end do
  end select
end subroutine gsstep_bwd_shift_dense

!!!!

subroutine gsstep_bwd_shift_csr(trans, n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_bwd_shift_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
        sigma, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_bwd_shift_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
        sigma, omega, b, incb, x, incx)
  end select
end subroutine gsstep_bwd_shift_csr

subroutine gsstep_bwd_shift_csc(trans, n, alpha, spA, colptr, rowind, nnz, &
  sigma, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_bwd_shift_csr_trans(n, alpha, spA, colptr, rowind, nnz, &
        sigma, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_bwd_shift_csr_notrans(n, alpha, spA, colptr, rowind, nnz, &
        sigma, omega, b, incb, x, incx)
  end select
end subroutine gsstep_bwd_shift_csc

!!! private

subroutine gsstep_bwd_shift_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  use sparse
  integer, parameter :: base = sparse_base_index
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(base:*), b(1:incb,base:base+n-1)
  integer, intent(in) :: rowptr(base:*), colind(base:*)
  double precision, intent(inout) :: x(1:incx,base:base+n-1)

  integer :: i, j, z
  double precision :: tmpd, tmpx, c

  do i = base+n-1, base, -1
    tmpd = 0.0d0
    c = 0.0d0
    tmpx = b(1,i) / alpha
    do z = rowptr(i+1)-1, rowptr(i), -1
      j = colind(z)
      if (i == j) then
        tmpd = spA(z)
!        tmpx = tmpx + sigma * x(1,i)
        call kahan_sum(tmpx, sigma * x(1,i), c)
      else
!        tmpx = tmpx - spA(z) * x(1,j)
        call kahan_sum(tmpx, - spA(z) * x(1,j), c)
      end if
    end do
    x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
  end do
end subroutine gsstep_bwd_shift_csr_notrans

subroutine gsstep_bwd_shift_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
  sigma, omega, b, incb, x, incx)
  use sparse
  integer, parameter :: base = sparse_base_index
  integer, intent(in) :: n, nnz, incb, incx
  double precision, intent(in) :: alpha, sigma, omega, spA(base:*), b(1:incb,base:base+n-1)
  integer, intent(in) :: rowptr(base:*), colind(base:*)
  double precision, intent(inout) :: x(1:incx,base:base+n-1)

  integer :: i, j, z
  double precision :: tmp(base:base+n-1)
  double precision :: c(base:base+n-1)
  integer :: diag(base:base+n-1)

  tmp = 0.0d0
  c = 0.0d0
  call daxpy(n, 1.0/alpha, b, incb, tmp, 1)
  do i = base+n-1, base, -1
    do z = rowptr(i+1)-1, rowptr(i), -1
      j = colind(z)
      if (i == j) then
        diag(i) = z
!        tmp(i) = tmp(i) + sigma * x(1,i)
        call kahan_sum(tmp(i), sigma * x(1,i), c(i))
        exit
      else
!        tmp(j) = tmp(j) - spA(z) * x(1,i)
        call kahan_sum(tmp(j), - spA(z) * x(1,i), c(j))
      end if
    end do
  end do

  do i = base+n-1, base, -1
    x(1,i) = omega * tmp(i) / spA(diag(i)) + (1.0d0 - omega) * x(1,i)
    do z = diag(i)-1, rowptr(i), -1
      j = colind(z)
!      tmp(j) = tmp(j) - spA(z) * x(1,i)
      call kahan_sum(tmp(j), - spA(z) * x(1,i), c(j))
    end do
  end do
end subroutine gsstep_bwd_shift_csr_trans

!!!!! matrix

subroutine gsstep_fwd_shift_dense_mat(trans, n, nrhs, alpha, A, lda, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, lda, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, A(1:lda,1:n), b(1:ldb,1:n)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_fwd_shift_dense(trans, n, alpha, A, lda, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_fwd_shift_dense_mat

subroutine gsstep_fwd_shift_csr_mat(trans, n, nrhs, alpha, spA, rowptr, colind, nnz, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:ldb,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_fwd_shift_csr(trans, n, alpha, spA, rowptr, colind, nnz, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_fwd_shift_csr_mat

subroutine gsstep_fwd_shift_csc_mat(trans, n, nrhs, alpha, spA, colptr, rowind, nnz, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:ldb,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_fwd_shift_csc(trans, n, alpha, spA, colptr, rowind, nnz, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_fwd_shift_csc_mat

!!!!! matrix

subroutine gsstep_bwd_shift_dense_mat(trans, n, nrhs, alpha, A, lda, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, lda, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, A(1:lda,1:n), b(1:ldb,1:n)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_bwd_shift_dense(trans, n, alpha, A, lda, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_bwd_shift_dense_mat

subroutine gsstep_bwd_shift_csr_mat(trans, n, nrhs, alpha, spA, rowptr, colind, nnz, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:ldb,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_bwd_shift_csr(trans, n, alpha, spA, rowptr, colind, nnz, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_bwd_shift_csr_mat

subroutine gsstep_bwd_shift_csc_mat(trans, n, nrhs, alpha, spA, colptr, rowind, nnz, sigma, omega, b, ldb, x, ldx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldb, ldx
  double precision, intent(in) :: alpha, sigma, omega, spA(1:nnz), b(1:ldb,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(inout) :: x(1:ldx,1:n)

  integer :: i

  do i = 1, nrhs
    call gsstep_bwd_shift_csc(trans, n, alpha, spA, colptr, rowind, nnz, sigma, omega, b(1,i), 1, x(1,i), 1)
  end do
end subroutine gsstep_bwd_shift_csc_mat

end module gsstep
