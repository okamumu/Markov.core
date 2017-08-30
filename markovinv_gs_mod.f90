
module markovinv_gs
  implicit none
  private
  public markovinv_gs_dense
  public markovinv_gs_csr
  public markovinv_gs_csc
!   public markovinv_gs_coo

contains

  ! Description: sojourn time for lossy CTMC
  !
  !        y = x * (-Q)^-1
  !
  !        or
  !
  !        y = (-Q)^-1 * x
  !
  ! Q: CTMC kernel

  subroutine markovinv_gs_dense(trans, n, nrhs, Q, ldq, x, ldx, &
    ystart, ldys, y, ldy, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, ldq, ldx, ldys, ldy, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:ldq,1:n), x(1:ldx,1:nrhs), ystart(1:ldys,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k, i
    double precision :: prevy(1:n,1:nrhs)

    y(1:n,1:nrhs) = ystart(1:n,1:nrhs)
    prevy(1:n,1:nrhs) = 0.0d0

    iter = 0
    do
      do i = 1, nrhs
        call dcopy(n, y(1,i), 1, prevy(1,i), 1)
        do k = 1, steps
          call gsstep_fwd_shift_dense(trans, n, -1.0d0, Q, ldq, 0.0d0, 1.0d0, x(1,i), 1, y(1,i), 1)
        end do
        call daxpy(n, -1.0d0, y(1,i), 1, prevy(1,i), 1)
      end do
      iter = iter + steps

      rerror = maxval(abs(prevy(1:n,1:nrhs)/y(1:n,1:nrhs)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovinv_gs_dense

  subroutine markovinv_gs_csr(trans, n, nrhs, Q, rowptr, colind, nnz, x, ldx, &
    ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, nnz, ldx, ldys, ldy, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:nnz), x(1:ldx,1:nrhs), ystart(1:ldys,1:nrhs)
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(out) :: y(1:ldy,1:nrhs), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k, i
    double precision :: prevy(1:n,1:nrhs)

    y(1:n,1:nrhs) = ystart(1:n,1:nrhs)
    prevy(1:n,1:nrhs) = 0.0d0

    iter = 0
    do
      do i = 1, nrhs
        call dcopy(n, y(1,i), 1, prevy(1,i), 1)
        do k = 1, steps
          call gsstep_fwd_shift_csr(trans, n, -1.0d0, Q, rowptr, colind, nnz, 0.0d0, 1.0d0, x(1,i), 1, y(1,i), 1)
        end do
        call daxpy(n, -1.0d0, y(1,i), 1, prevy(1,i), 1)
      end do
      iter = iter + steps

      rerror = maxval(abs(prevy(1:n,1:nrhs)/y(1:n,1:nrhs)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovinv_gs_csr

  subroutine markovinv_gs_csc(trans, n, nrhs, Q, colptr, rowind, nnz, x, ldx, &
    ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, nnz, ldx, ldys, ldy, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:nnz), x(1:ldx,1:nrhs), ystart(1:ldys,1:nrhs)
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(out) :: y(1:ldy,1:nrhs), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k, i
    double precision :: prevy(1:n,1:nrhs)

    y(1:n,1:nrhs) = ystart(1:n,1:nrhs)
    prevy(1:n,1:nrhs) = 0.0d0

    iter = 0
    do
      do i = 1, nrhs
        call dcopy(n, y(1,i), 1, prevy(1,i), 1)
        do k = 1, steps
          call gsstep_fwd_shift_csc(trans, n, -1.0d0, Q, colptr, rowind, nnz, 0.0d0, 1.0d0, x(1,i), 1, y(1,i), 1)
        end do
        call daxpy(n, -1.0d0, y(1,i), 1, prevy(1,i), 1)
      end do
      iter = iter + steps

      rerror = maxval(abs(prevy(1:n,1:nrhs)/y(1:n,1:nrhs)))

      ! callback
      call callback(iter, rerror)

      if (rerror < rtol) then
        info = 0
        exit
      end if

      if (iter >= maxiter) then
        info = -1
        exit
      end if
    end do
  end subroutine markovinv_gs_csc

end module markovinv_gs
