
module markovstsen_gs
  implicit none
  private
  public markovstsen_gs_dense
  public markovstsen_gs_csr
  public markovstsen_gs_csc
!   public markovstsen_gs_coo

contains

  ! Description: sensitivity vector for stationary vector
  !
  !        s * Q + pis * dQ = 0,  s*1 = 0
  !
  ! Q: CTMC kernel
  ! dQ: the first derivative of CTMC kernel with respect to a parameter
  ! pis: stationary vector such as pis * Q = 0, pis*1 = 1

  subroutine markovstsen_gs_dense(n, Q, ldq, dQ, lddq, pis, incp, &
    sstart, incss, s, incs, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use gsstep
    integer, intent(in) :: n, ldq, lddq, incp, incss, incs, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n), pis(1:incp,1:n), sstart(1:incss,1:n)
    double precision, intent(out) :: s(1:incs,1:n), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: tmp, prevs(1:n), b(1:n)

    s(1,1:n) = sstart(1,1:n)
    prevs(1:n) = 0.0d0
    call dgemv('T', n, n, 1.0d0, dQ, lddq, pis, incp, 0.0d0, b, 1)

    iter = 0
    do
      call dcopy(n, s, incs, prevs, 1)
      do k = 1, steps
        call gsstep_fwd_shift_dense('T', n, -1.0d0, Q, ldq, 0.0d0, 1.0d0, b, 1, s, incs)
        tmp = sum(s(1,1:n))
        call daxpy(n, -tmp, pis, incp, s, incs)
      end do
      call daxpy(n, -1.0d0, s, incs, prevs, 1)
      iter = iter + steps

      rerror = maxval(abs(prevs(1:n)/s(1,1:n)))

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
  end subroutine markovstsen_gs_dense

  subroutine markovstsen_gs_csr(n, Q, rowptr0, colind0, nnz0, &
    dQ, rowptr1, colind1, nnz1, pis, incp, &
    sstart, incss, s, incs, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:nnz0), dQ(1:nnz1), pis(1:incp,1:n), sstart(1:incss,1:n)
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(out) :: s(1:incs,1:n), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: tmp, prevs(1:n), b(1:n)

    s(1,1:n) = sstart(1,1:n)
    prevs(1:n) = 0.0d0
    call dcsrmv('T', n, n, 1.0d0, dQ, rowptr1, colind1, nnz1, pis, incp, 0.0d0, b, 1)

    iter = 0
    do
      call dcopy(n, s, incs, prevs, 1)
      do k = 1, steps
        call gsstep_fwd_shift_csr('T', n, -1.0d0, Q, rowptr0, colind0, nnz0, 0.0d0, 1.0d0, b, 1, s, incs)
        tmp = sum(s(1,1:n))
        call daxpy(n, -tmp, pis, incp, s, incs)
      end do
      call daxpy(n, -1.0d0, s, incs, prevs, 1)
      iter = iter + steps

      rerror = maxval(abs(prevs(1:n)/s(1,1:n)))

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
  end subroutine markovstsen_gs_csr

  subroutine markovstsen_gs_csc(n, Q, colptr0, rowind0, nnz0, &
    dQ, colptr1, rowind1, nnz1, pis, incp, &
    sstart, incss, s, incs, &
    maxiter, rtol, steps, iter, rerror, info, callback)
    use spblas
    use gsstep
    integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
    double precision :: rtol
    double precision, intent(in) :: Q(1:nnz0), dQ(1:nnz1), pis(1:incp,1:n), sstart(1:incss,1:n)
    integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
    integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
    double precision, intent(out) :: s(1:incs,1:n), rerror
    integer, intent(out) :: iter, info

    interface
      subroutine callback(iter, rerror)
        integer, intent(in) :: iter
        double precision, intent(in) :: rerror
      end subroutine callback
    end interface

    integer :: k
    double precision :: tmp, prevs(1:n), b(1:n)

    s(1,1:n) = sstart(1,1:n)
    prevs(1:n) = 0.0d0
    call dcscmv('T', n, n, 1.0d0, dQ, colptr1, rowind1, nnz1, pis, incp, 0.0d0, b, 1)

    iter = 0
    do
      call dcopy(n, s, incs, prevs, 1)
      do k = 1, steps
        call gsstep_fwd_shift_csc('T', n, -1.0d0, Q, colptr0, rowind0, nnz0, 0.0d0, 1.0d0, b, 1, s, incs)
        tmp = sum(s(1,1:n))
        call daxpy(n, -tmp, pis, incp, s, incs)
      end do
      call daxpy(n, -1.0d0, s, incs, prevs, 1)
      iter = iter + steps

      rerror = maxval(abs(prevs(1:n)/s(1,1:n)))

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
  end subroutine markovstsen_gs_csc

end module markovstsen_gs
