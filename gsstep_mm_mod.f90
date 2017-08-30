!
! gsstep shift
!

module gsstep_mm
  implicit none
  public gsstep_fwd_mm_dense, gsstep_fwd_mm_csr, gsstep_fwd_mm_csc
  public gsstep_bwd_mm_dense, gsstep_bwd_mm_csr, gsstep_bwd_mm_csc
  private gsstep_fwd_mm_csr_notrans, gsstep_fwd_mm_csr_trans
  private gsstep_bwd_mm_csr_notrans, gsstep_bwd_mm_csr_trans

contains

! Gauss-Seidel step for solving the following linear equation
!
!         alpha * trans(A + f * g) * x = b
!
!         The step computes
!
!         x := (D + L)^(-1) (b/alpha - U) * x
!
!         A: square matrix
!         x: vector (in; initial vector for the step, out; updated vector)
!         b: constant vector
!         f: constant (column) vector
!         g: constant (row) vector

  subroutine gsstep_fwd_mm_dense(trans, n, alpha, A, lda, f, incf, g, incg, omega, b, incb, x, incx)
    character, intent(in) :: trans
    integer, intent(in) :: n, lda, incb, incx, incf, incg
    double precision, intent(in) :: alpha, omega, A(1:lda,1:n), b(1:incb,1:n)
    double precision, intent(inout) :: x(1:incx,1:n)
    double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

    integer :: i, j
    double precision :: tmpx, tmpd

    select case (trans)
      case ('N', 'n')
        do i = 1, n
          tmpd = 0.0d0
          tmpx = b(1,i) / alpha
          do j = 1, n
            if (i == j) then
              tmpd = A(i,j) + f(1,i) * g(1,j)
            else
              tmpx = tmpx - (A(i,j) + f(1,i) * g(1,j)) * x(1,j)
            end if
          end do
          x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
        end do
      case ('T', 't')
        do j = 1, n
          tmpd = 0.0d0
          tmpx = b(1,j) / alpha
          do i = 1, n
            if (i == j) then
              tmpd = A(i,j) + f(1,i) * g(1,j)
            else
              tmpx = tmpx - (A(i,j) + f(1,i) * g(1,j)) * x(1,i)
            end if
          end do
          x(1,j) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,j)
        end do
    end select
  end subroutine gsstep_fwd_mm_dense

  !!!!

  subroutine gsstep_fwd_mm_csr(trans, n, alpha, spA, rowptr, colind, nnz, &
    f, incf, g, incg, omega, b, incb, x, incx)
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, incb, incx, incf, incg
    double precision, intent(in) :: alpha, omega, spA(1:nnz), b(1:incb,1:n)
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(inout) :: x(1:incx,1:n)
    double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

    select case (trans)
      case ('N', 'n')
        call gsstep_fwd_mm_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
          f, incf, g, incg, omega, b, incb, x, incx)
      case ('T', 't')
        call gsstep_fwd_mm_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
          f, incf, g, incg, omega, b, incb, x, incx)
    end select
  end subroutine gsstep_fwd_mm_csr

  subroutine gsstep_fwd_mm_csc(trans, n, alpha, spA, colptr, rowind, nnz, &
    f, incf, g, incg, omega, b, incb, x, incx)
    character, intent(in) :: trans
    integer, intent(in) :: n, nnz, incb, incx, incf, incg
    double precision, intent(in) :: alpha, omega, spA(1:nnz), b(1:incb,1:n)
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(inout) :: x(1:incx,1:n)
    double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

    select case (trans)
      case ('N', 'n')
        call gsstep_fwd_mm_csr_trans(n, alpha, spA, colptr, rowind, nnz, &
          g, incg, f, incf, omega, b, incb, x, incx)
      case ('T', 't')
        call gsstep_fwd_mm_csr_notrans(n, alpha, spA, colptr, rowind, nnz, &
          g, incg, f, incf, omega, b, incb, x, incx)
    end select
  end subroutine gsstep_fwd_mm_csc

  subroutine gsstep_fwd_mm_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
      f, incf, g, incg, omega, b, incb, x, incx)
      use sparse
      integer, parameter :: base = sparse_base_index
      integer, intent(in) :: n, nnz, incb, incx, incf, incg
      double precision, intent(in) :: alpha, omega, spA(base:*), b(1:incb,base:base+n-1)
      integer, intent(in) :: rowptr(base:*), colind(base:*)
      double precision, intent(inout) :: x(1:incx,base:base+n-1)
      double precision, intent(in) :: f(1:incf,base:base+n-1), g(1:incg,base:base+n-1)

      integer :: i, j, z
      double precision :: tmpd, tmpx
      double precision :: sumxg
      double precision :: xg(base:base+n-1)

      xg(base:base+n-1) = x(1,base:base+n-1) * g(1,base:base+n-1)
      sumxg = sum(xg(base+1:base+n-1))
      do i = base, base+n-1
        tmpd = 0.0d0
        tmpx = b(1,i) / alpha
        if (i /= base) then
          sumxg = sumxg - xg(i) + x(1,i-1) * g(1,i-1)
        else
          sumxg = sumxg - xg(i)
        end if
        tmpx = tmpx - f(1,i) * sumxg
        do z = rowptr(i), rowptr(i+1)-1
          j = colind(z)
          if (i == j) then
            tmpd = spA(z) + f(1,i) * g(1,j)
          else
            tmpx = tmpx - spA(z) * x(1,j)
          end if
        end do
        x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
      end do
    end subroutine gsstep_fwd_mm_csr_notrans

    subroutine gsstep_fwd_mm_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
      f, incf, g, incg, omega, b, incb, x, incx)
      use sparse
      integer, parameter :: base = sparse_base_index
      integer, intent(in) :: n, nnz, incb, incx, incf, incg
      double precision, intent(in) :: alpha, omega, spA(base:*), b(1:incb,base:base+n-1)
      integer, intent(in) :: rowptr(base:*), colind(base:*)
      double precision, intent(inout) :: x(1:incx,base:base+n-1)
      double precision, intent(in) :: f(1:incf,base:base+n-1), g(1:incg,base:base+n-1)

      integer :: i, j, z
      double precision :: tmp(base:base+n-1)
      integer :: diag(base:base+n-1)
      double precision :: sumxg
      double precision :: xg(base:base+n-1)

      tmp = 0.0d0
      call daxpy(n, 1.0/alpha, b, incb, tmp, 1)
      do i = base, base+n-1
        do z = rowptr(i), rowptr(i+1)-1
          j = colind(z)
          if (i == j) then
            diag(i) = z
            exit
          else
            tmp(j) = tmp(j) - spA(z) * x(1,i)
          end if
        end do
      end do

      xg(base:base+n-1) = x(1,base:base+n-1) * f(1,base:base+n-1)
      sumxg = sum(xg(base+1:base+n-1))

      do i = base, base+n-1
        if (i /= base) then
          sumxg = sumxg - xg(i) + x(1,i-1) * f(1,i-1)
        else
          sumxg = sumxg - xg(i)
        end if
        tmp(i) = tmp(i) - g(1,i) * sumxg
        x(1,i) = omega * tmp(i) / (spA(diag(i)) + f(1,i) * g(1,i)) + (1.0d0 - omega) * x(1,i)
        do z = diag(i)+1, rowptr(i+1)-1
          j = colind(z)
          tmp(j) = tmp(j) - spA(z) * x(1,i)
        end do
      end do
    end subroutine gsstep_fwd_mm_csr_trans

!!!! backward

subroutine gsstep_bwd_mm_dense(trans, n, alpha, A, lda, f, incf, g, incg, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, incb, incx, incf, incg
  double precision, intent(in) :: alpha, omega, A(1:lda,1:n), b(1:incb,1:n)
  double precision, intent(inout) :: x(1:incx,1:n)
  double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

  integer :: i, j
  double precision :: tmpx, tmpd

  select case (trans)
    case ('N', 'n')
      do i = n, 1, -1
        tmpd = 0.0d0
        tmpx = b(1,i) / alpha
        do j = n, 1, -1
          if (i == j) then
            tmpd = A(i,j) + f(1,i) * g(1,j)
          else
            tmpx = tmpx - (A(i,j) + f(1,i) * g(1,j)) * x(1,j)
          end if
        end do
        x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
      end do
    case ('T', 't')
      do j = n, 1, -1
        tmpd = 0.0d0
        tmpx = b(1,j) / alpha
        do i = n, 1, -1
          if (i == j) then
            tmpd = A(i,j) + f(1,i) * g(1,j)
          else
            tmpx = tmpx - (A(i,j) + f(1,i) * g(1,j)) * x(1,i)
          end if
        end do
        x(1,j) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,j)
      end do
  end select
end subroutine gsstep_bwd_mm_dense

!!!!

subroutine gsstep_bwd_mm_csr(trans, n, alpha, spA, rowptr, colind, nnz, &
  f, incf, g, incg, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx, incf, incg
  double precision, intent(in) :: alpha, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)
  double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_bwd_mm_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
        f, incf, g, incg, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_bwd_mm_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
        f, incf, g, incg, omega, b, incb, x, incx)
  end select
end subroutine gsstep_bwd_mm_csr

subroutine gsstep_bwd_mm_csc(trans, n, alpha, spA, colptr, rowind, nnz, &
  f, incf, g, incg, omega, b, incb, x, incx)
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incb, incx, incf, incg
  double precision, intent(in) :: alpha, omega, spA(1:nnz), b(1:incb,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(inout) :: x(1:incx,1:n)
  double precision, intent(in) :: f(1:incf,1:n), g(1:incg,1:n)

  select case (trans)
    case ('N', 'n')
      call gsstep_bwd_mm_csr_trans(n, alpha, spA, colptr, rowind, nnz, &
        g, incg, f, incf, omega, b, incb, x, incx)
    case ('T', 't')
      call gsstep_bwd_mm_csr_notrans(n, alpha, spA, colptr, rowind, nnz, &
        g, incg, f, incf, omega, b, incb, x, incx)
  end select
end subroutine gsstep_bwd_mm_csc

subroutine gsstep_bwd_mm_csr_notrans(n, alpha, spA, rowptr, colind, nnz, &
    f, incf, g, incg, omega, b, incb, x, incx)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz, incb, incx, incf, incg
    double precision, intent(in) :: alpha, omega, spA(base:*), b(1:incb,base:base+n-1)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    double precision, intent(inout) :: x(1:incx,base:base+n-1)
    double precision, intent(in) :: f(1:incf,base:base+n-1), g(1:incg,base:base+n-1)

    integer :: i, j, z
    double precision :: tmpd, tmpx
    double precision :: sumxg
    double precision :: xg(base:base+n-1)

    xg(base:base+n-1) = x(1,base:base+n-1) * g(1,base:base+n-1)
    sumxg = sum(xg(base+1:base+n-1))
    do i = base+n-1, base, -1
      tmpd = 0.0d0
      tmpx = b(1,i) / alpha
      if (i /= base) then
        sumxg = sumxg - xg(i) + x(1,i+1) * g(1,i+1)
      else
        sumxg = sumxg - xg(i)
      end if
      tmpx = tmpx - f(1,i) * sumxg
      do z = rowptr(i+1)-1, rowptr(i), -1
        j = colind(z)
        if (i == j) then
          tmpd = spA(z) + f(1,i) * g(1,j)
        else
          tmpx = tmpx - spA(z) * x(1,j)
        end if
      end do
      x(1,i) = omega * tmpx / tmpd + (1.0d0 - omega) * x(1,i)
    end do
  end subroutine gsstep_bwd_mm_csr_notrans

  subroutine gsstep_bwd_mm_csr_trans(n, alpha, spA, rowptr, colind, nnz, &
    f, incf, g, incg, omega, b, incb, x, incx)
    use sparse
    integer, parameter :: base = sparse_base_index
    integer, intent(in) :: n, nnz, incb, incx, incf, incg
    double precision, intent(in) :: alpha, omega, spA(base:*), b(1:incb,base:base+n-1)
    integer, intent(in) :: rowptr(base:*), colind(base:*)
    double precision, intent(inout) :: x(1:incx,base:base+n-1)
    double precision, intent(in) :: f(1:incf,base:base+n-1), g(1:incg,base:base+n-1)

    integer :: i, j, z
    double precision :: tmp(base:base+n-1)
    integer :: diag(base:base+n-1)
    double precision :: sumxg
    double precision :: xg(base:base+n-1)

    tmp = 0.0d0
    call daxpy(n, 1.0/alpha, b, incb, tmp, 1)
    do i = base+n-1, base, -1
      do z = rowptr(i+1)-1, rowptr(i), -1
        j = colind(z)
        if (i == j) then
          diag(i) = z
          exit
        else
          tmp(j) = tmp(j) - spA(z) * x(1,i)
        end if
      end do
    end do

    xg(base:base+n-1) = x(1,base:base+n-1) * f(1,base:base+n-1)
    sumxg = sum(xg(base+1:base+n-1))

    do i = base+n-1, base
      if (i /= base) then
        sumxg = sumxg - xg(i) + x(1,i+1) * f(1,i+1)
      else
        sumxg = sumxg - xg(i)
      end if
      tmp(i) = tmp(i) - g(1,i) * sumxg
      x(1,i) = omega * tmp(i) / (spA(diag(i)) + f(1,i) * g(1,i)) + (1.0d0 - omega) * x(1,i)
      do z = diag(i)-1, rowptr(i), -1
        j = colind(z)
        tmp(j) = tmp(j) - spA(z) * x(1,i)
      end do
    end do
  end subroutine gsstep_bwd_mm_csr_trans
end module gsstep_mm
