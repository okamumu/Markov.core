!
! f90 wrapper
!

subroutine f90_markovinv_dgesv_dense(trans, n, nrhs, Q, ldq, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, ldq, ldx, ldy
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(in) :: x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  ! to avoid an error caused by a bug in gfortran?
  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_dense(trans, n, nrhs, Q, ldq, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_dense

subroutine f90_markovinv_dgesv_csr(trans, n, nrhs, spQ, rowptr, colind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(in) :: x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_csr(trans, n, nrhs, spQ, rowptr, colind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_csr

subroutine f90_markovinv_dgesv_csc(trans, n, nrhs, spQ, colptr, rowind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_csc(trans, n, nrhs, spQ, colptr, rowind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_csc

subroutine f90_markovinv_dgesv_coo(trans, n, nrhs, spQ, rowind, colind, nnz, x, ldx, y, ldy, info)
  use markovinv_dgesv
  character, intent(in) :: trans
  integer, intent(in) :: n, nrhs, nnz, ldx, ldy
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
  double precision, intent(out) :: y(1:ldy,1:nrhs)
  integer, intent(out) :: info

  integer :: tmp
  tmp = ldy
  call markovinv_dgesv_coo(trans, n, nrhs, spQ, rowind, colind, nnz, x, ldx, y, tmp, info)
end subroutine f90_markovinv_dgesv_coo

!!!!

subroutine f90_markovinv_gs_dense(trans, n, nrhs, Q, ldq, x, ldx, &
  ystart, ldys, y, ldy, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovinv_gs
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

  call markovinv_gs_dense(trans, n, nrhs, Q, ldq, x, ldx, &
    ystart, ldys, y, ldy, &
    maxiter, rtol, steps, iter, rerror, info, callback)

end subroutine f90_markovinv_gs_dense

subroutine f90_markovinv_gs_csr(trans, n, nrhs, Q, rowptr, colind, nnz, x, ldx, &
  ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovinv_gs
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

  call markovinv_gs_csr(trans, n, nrhs, Q, rowptr, colind, nnz, x, ldx, &
    ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)

end subroutine f90_markovinv_gs_csr

subroutine f90_markovinv_gs_csc(trans, n, nrhs, Q, colptr, rowind, nnz, x, ldx, &
  ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovinv_gs
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

  call markovinv_gs_csc(trans, n, nrhs, Q, colptr, rowind, nnz, x, ldx, &
    ystart, ldys, y, ldy, maxiter, rtol, steps, iter, rerror, info, callback)

end subroutine f90_markovinv_gs_csc

!!!

subroutine f90_markovqst2_gs_dense(n, Q, ldq, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, ldq, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_dense(n, Q, ldq, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_dense

subroutine f90_markovqst2_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, nnz, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_csr

subroutine f90_markovqst2_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
  use markovqst2_gs
  integer, intent(in) :: n, nnz, incxi, incxs, incx, incys, incy, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n)
  double precision, intent(in) :: xstart(1:incxs,1:n), ystart(1:incys,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: gam, rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovqst2_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
  xstart, incxs, x, incx, ystart, incys, y, incy, &
  gam, maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovqst2_gs_csc

!!!

! subroutine f90_markovqst_gs_dense(n, Q, ldq, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
!   use markovqst_gs
!   integer, intent(in) :: n, ldq, incxi, incxs, incx, maxiter, steps
!   double precision, intent(in) :: rtol
!   double precision, intent(in) :: Q(1:ldq,1:n), xi(1:incxi,1:n), xstart(1:incxs,1:n)
!   double precision, intent(out) :: x(1:incx,1:n), gam, rerror
!   integer, intent(out) :: iter, info

!   interface
!     subroutine callback(iter, rerror)
!       integer, intent(in) :: iter
!       double precision, intent(in) :: rerror
!     end subroutine callback
!   end interface

!   call markovqst_gs_dense(n, Q, ldq, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
! end subroutine f90_markovqst_gs_dense

! subroutine f90_markovqst_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
!   use markovqst_gs
!   integer, intent(in) :: n, nnz, incxi, incxs, incx, maxiter, steps
!   double precision, intent(in) :: rtol
!   double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n), xstart(1:incxs,1:n)
!   integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
!   double precision, intent(out) :: x(1:incx,1:n), gam, rerror
!   integer, intent(out) :: iter, info

!   interface
!     subroutine callback(iter, rerror)
!       integer, intent(in) :: iter
!       double precision, intent(in) :: rerror
!     end subroutine callback
!   end interface

!   call markovqst_gs_csr(n, Q, rowptr, colind, nnz, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
! end subroutine f90_markovqst_gs_csr

! subroutine f90_markovqst_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
!   use markovqst_gs
!   integer, intent(in) :: n, nnz, incxi, incxs, incx, maxiter, steps
!   double precision, intent(in) :: rtol
!   double precision, intent(in) :: Q(1:nnz), xi(1:incxi,1:n), xstart(1:incxs,1:n)
!   integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
!   double precision, intent(out) :: x(1:incx,1:n), gam, rerror
!   integer, intent(out) :: iter, info

!   interface
!     subroutine callback(iter, rerror)
!       integer, intent(in) :: iter
!       double precision, intent(in) :: rerror
!     end subroutine callback
!   end interface

!   call markovqst_gs_csc(n, Q, colptr, rowind, nnz, xi, incxi, &
!   xstart, incxs, x, incx, &
!   gam, maxiter, rtol, steps, iter, rerror, info, callback)
! end subroutine f90_markovqst_gs_csc

!!!

! subroutine f90_markovst_dgesv_dense(n, Q, ldq, x, incx, info)
!   use markovst_dgesv
!   integer, intent(in) :: n, ldq, incx
!   double precision, intent(in) :: Q(1:ldq,1:n)
!   double precision, intent(out) :: x(1:incx,1:n)
!   integer, intent(out) :: info

!   call markovst_dgesv_dense(n, Q, ldq, x, incx, info)
! end subroutine f90_markovst_dgesv_dense

! subroutine f90_markovst_dgesv_csr(n, spQ, rowptr, colind, nnz, x, incx, info)
!   use markovst_dgesv
!   integer, intent(in) :: n, nnz, incx
!   integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
!   double precision, intent(in) :: spQ(1:nnz)
!   double precision, intent(out) :: x(1:incx,1:n)
!   integer, intent(out) :: info

!   call markovst_dgesv_csr(n, spQ, rowptr, colind, nnz, x, incx, info)
! end subroutine f90_markovst_dgesv_csr

! subroutine f90_markovst_dgesv_csc(n, spQ, colptr, rowind, nnz, x, incx, info)
!   use markovst_dgesv
!   integer, intent(in) :: n, nnz, incx
!   integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
!   double precision, intent(in) :: spQ(1:nnz)
!   double precision, intent(out) :: x(1:incx,1:n)
!   integer, intent(out) :: info

!   call markovst_dgesv_csc(n, spQ, colptr, rowind, nnz, x, incx, info)
! end subroutine f90_markovst_dgesv_csc

! subroutine f90_markovst_dgesv_coo(n, spQ, rowind, colind, nnz, x, incx, info)
!   use markovst_dgesv
!   integer, intent(in) :: n, nnz, incx
!   integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
!   double precision, intent(in) :: spQ(1:nnz)
!   double precision, intent(out) :: x(1:incx,1:n)
!   integer, intent(out) :: info

!   call markovst_dgesv_coo(n, spQ, rowind, colind, nnz, x, incx, info)
! end subroutine f90_markovst_dgesv_coo

!!!

subroutine f90_markovst_gs_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_dense

subroutine f90_markovst_gs_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_csr

subroutine f90_markovst_gs_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovst_gs
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovst_gs_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovst_gs_csc

!!!

subroutine f90_markovst_sor_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, ldq, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), xstart(1:incxs,1:n)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_dense(n, Q, ldq, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_dense

subroutine f90_markovst_sor_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_csr(n, Q, rowptr, colind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_csr

subroutine f90_markovst_sor_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
  use markovst_sor
  integer, intent(in) :: n, nnz, incxs, incx, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:nnz), xstart(1:incxs,1:n)
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n), rerror
  double precision, intent(inout) :: omega
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, steps, rerror, omega, n, prevx, x, incx)
      integer, intent(in) :: iter, n, incx
      integer, intent(inout) :: steps
      double precision, intent(in) :: rerror, omega
      double precision, intent(in) :: prevx(1:n), x(1:incx,1:n)
    end subroutine callback
  end interface

  interface
    subroutine update_omega(iter, n, x1, x2, omega)
      integer, intent(in) :: iter, n
      double precision, intent(in) :: x1(1:n), x2(1:n)
      double precision, intent(inout) :: omega
    end subroutine update_omega
  end interface

  call markovst_sor_csc(n, Q, colptr, rowind, nnz, xstart, incxs, x, incx, &
  maxiter, rtol, steps, iter, rerror, omega, info, callback, update_omega)
end subroutine f90_markovst_sor_csc

!!!

subroutine f90_markovst_gth_dense(n, Q, ldq, x, incx)
  use markovst_gth
  integer, intent(in) :: n, ldq, incx
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_dense(n, Q, ldq, x, incx)
end subroutine f90_markovst_gth_dense

subroutine f90_markovst_gth_csr(n, spQ, rowptr, colind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_csr(n, spQ, rowptr, colind, nnz, x, incx)
end subroutine f90_markovst_gth_csr

subroutine f90_markovst_gth_csc(n, spQ, colptr, rowind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_csc(n, spQ, colptr, rowind, nnz, x, incx)
end subroutine f90_markovst_gth_csc

subroutine f90_markovst_gth_coo(n, spQ, rowind, colind, nnz, x, incx)
  use markovst_gth
  integer, intent(in) :: n, nnz, incx
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: x(1:incx,1:n)

  call markovst_gth_coo(n, spQ, rowind, colind, nnz, x, incx)
end subroutine f90_markovst_gth_coo

!!!

subroutine f90_markovstsen_dgesv_dense(n, Q, ldq, dQ, lddq, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, ldq, lddq, incp, incs
  double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_dense(n, Q, ldq, dQ, lddq, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_dense

subroutine f90_markovstsen_dgesv_csr(n, spQ, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
  integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_csr(n, spQ, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_csr

subroutine f90_markovstsen_dgesv_csc(n, spQ, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
  integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_csc(n, spQ, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_csc

subroutine f90_markovstsen_dgesv_coo(n, spQ, rowind0, colind0, nnz0, &
  dQ, rowind1, colind1, nnz1, pis, incp, s, incs, info)
  use markovstsen_dgesv
  integer, intent(in) :: n, nnz0, nnz1, incp, incs
  integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
  integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
  double precision, intent(in) :: spQ(1:nnz0)
  double precision, intent(in) :: dQ(1:nnz1)
  double precision, intent(in) :: pis(1:incp,1:n)
  double precision, intent(out) :: s(1:incs,1:n)
  integer, intent(out) :: info

  call markovstsen_dgesv_coo(n, spQ, rowind0, colind0, nnz0, &
  dQ, rowind1, colind1, nnz1, pis, incp, s, incs, info)
end subroutine f90_markovstsen_dgesv_coo

!!!

subroutine f90_markovstsen_gs_dense(n, Q, ldq, dQ, lddq, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, ldq, lddq, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
  double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n), pis(1:incp,1:n), sstart(1:incss,1:n)
  double precision, intent(out) :: s(1:incs,1:n), rerror
  integer, intent(out) :: iter, info

  interface
    subroutine callback(iter, rerror)
      integer, intent(in) :: iter
      double precision, intent(in) :: rerror
    end subroutine callback
  end interface

  call markovstsen_gs_dense(n, Q, ldq, dQ, lddq, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_dense

subroutine f90_markovstsen_gs_csr(n, Q, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
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

  call markovstsen_gs_csr(n, Q, rowptr0, colind0, nnz0, &
  dQ, rowptr1, colind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_csr

subroutine f90_markovstsen_gs_csc(n, Q, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
  use markovstsen_gs
  integer, intent(in) :: n, nnz0, nnz1, incp, incss, incs, maxiter, steps
  double precision, intent(in) :: rtol
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

  call markovstsen_gs_csc(n, Q, colptr0, rowind0, nnz0, &
  dQ, colptr1, rowind1, nnz1, pis, incp, &
  sstart, incss, s, incs, &
  maxiter, rtol, steps, iter, rerror, info, callback)
end subroutine f90_markovstsen_gs_csc

!!!

subroutine f90_mexp_pade_dense(trans, n, alpha, MA, lda, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, lde
  double precision, intent(in) :: alpha, MA(1:lda,1:n), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_dense(trans, n, alpha, MA, lda, ME, lde, eps)
end subroutine f90_mexp_pade_dense

subroutine f90_mexp_pade_csr(trans, n, alpha, spMA, rowptr, colind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_csr(trans, n, alpha, spMA, rowptr, colind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_csr

subroutine f90_mexp_pade_csc(trans, n, alpha, spMA, colptr, rowind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_csc(trans, n, alpha, spMA, colptr, rowind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_csc

subroutine f90_mexp_pade_coo(trans, n, alpha, spMA, rowind, colind, nnz, ME, lde, eps)
  use mexp_pade
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: alpha, spMA(1:*), eps
  double precision, intent(out) :: ME(1:lde,1:n)

  call mexp_pade_coo(trans, n, alpha, spMA, rowind, colind, nnz, ME, lde, eps)
end subroutine f90_mexp_pade_coo

!!!

subroutine f90_mexp_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, incx, incy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), x(1:*), atol
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_dense_vec

subroutine f90_mexp_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_csr_vec

subroutine f90_mexp_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_csc_vec

subroutine f90_mexp_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n)

  call mexp_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, atol)
end subroutine f90_mexp_unif_coo_vec

subroutine f90_mexp_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, m, ldx, ldy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_dense_mat

subroutine f90_mexp_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_csr_mat

subroutine f90_mexp_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_csc_mat

subroutine f90_mexp_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
  use mexp_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m)

  call mexp_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, atol)
end subroutine f90_mexp_unif_coo_mat

!!!

subroutine f90_mexpconv_unif_dense_vec(transQ, transH, n, P, ldp, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, H, ldh, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, ldp, left, right, incx, incy, incz, ldh
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), H(1:ldh,1:n)

  call mexpconv_unif_dense_vec(transQ, transH, n, P, ldp, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, H, ldh, atol)
end subroutine f90_mexpconv_unif_dense_vec

subroutine f90_mexpconv_unif_csr_vec(transQ, transH, n, spP, rowptr, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_csr_vec(transQ, transH, n, spP, rowptr, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_csr_vec

subroutine f90_mexpconv_unif_csc_vec(transQ, transH, n, spP, colptr, rowind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_csc_vec(transQ, transH, n, spP, colptr, rowind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_csc_vec

subroutine f90_mexpconv_unif_coo_vec(transQ, transH, n, spP, rowind, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
  use mexpconv_unif
  character, intent(in) :: transQ, transH
  integer, intent(in) :: n, nnz, left, right, incx, incy, incz
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: x(1:incx,1:n), y(1:incy,1:n)
  double precision, intent(out) :: z(1:incz,1:n), spH(1:*)

  call mexpconv_unif_coo_vec(transQ, transH, n, spP, rowind, colind, nnz, &
  qv, left, right, poi, weight, x, incx, y, incy, z, incz, spH, atol)
end subroutine f90_mexpconv_unif_coo_vec

!!!

subroutine f90_mexpint_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, incx, incy, inccy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), x(1:*), atol
  double precision, intent(out) :: y(1:incy,1:n), cy(1:*)

  call mexpint_unif_dense_vec(trans, n, P, ldp, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
end subroutine f90_mexpint_unif_dense_vec

subroutine f90_mexpint_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy, inccy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n), cy(1:*)

  call mexpint_unif_csr_vec(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
end subroutine f90_mexpint_unif_csr_vec

subroutine f90_mexpint_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy, inccy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(out) :: y(1:incy,1:n), cy(1:*)

  call mexpint_unif_csc_vec(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
end subroutine f90_mexpint_unif_csc_vec

subroutine f90_mexpint_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, incx, incy, inccy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), x(1:*), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(out) :: y(1:incy,1:n), cy(1:*)

  call mexpint_unif_coo_vec(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, x, incx, y, incy, cy, inccy, atol)
end subroutine f90_mexpint_unif_coo_vec

subroutine f90_mexpint_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, ldp, left, right, m, ldx, ldy, ldcy
  double precision, intent(in) :: P(1:ldp,1:n), qv, weight, poi(left:right), atol
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m), cy(1:ldcy,1:m)

  call mexpint_unif_dense_mat(trans, n, P, ldp, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
end subroutine f90_mexpint_unif_dense_mat

subroutine f90_mexpint_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy, ldcy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m), cy(1:ldcy,1:m)

  call mexpint_unif_csr_mat(trans, n, spP, rowptr, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
end subroutine f90_mexpint_unif_csr_mat

subroutine f90_mexpint_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy, ldcy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m), cy(1:ldcy,1:m)

  call mexpint_unif_csc_mat(trans, n, spP, colptr, rowind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
end subroutine f90_mexpint_unif_csc_mat

subroutine f90_mexpint_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
  use mexpint_unif
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, left, right, m, ldx, ldy, ldcy
  double precision, intent(in) :: spP(1:*), qv, weight, poi(left:right), atol
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: x(1:ldx,1:m)
  double precision, intent(out) :: y(1:ldy,1:m), cy(1:ldcy,1:m)

  call mexpint_unif_coo_mat(trans, n, spP, rowind, colind, nnz, qv, &
  left, right, poi, weight, m, x, ldx, y, ldy, cy, ldcy, atol)
end subroutine f90_mexpint_unif_coo_mat

!!!

subroutine f90_poisson_rightbound(lambda, eps, right)
  use poisson
  double precision, intent(in) :: lambda, eps
  integer, intent(out) :: right

  right = poisson_rightbound(lambda, eps)
end subroutine f90_poisson_rightbound

subroutine f90_poisson_prob(lambda, left, right, prob, weight)
  use poisson
  integer, intent(in) :: left, right
  double precision, intent(in) :: lambda
  double precision, intent(out) :: prob(left:right), weight

  call poisson_prob(lambda, left, right, prob, weight)
end subroutine f90_poisson_prob

!!!

subroutine f90_unif_dense(n, Q, ldq, P, ldp, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, ldq, ldp
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: Q(1:ldq,1:n)
  double precision, intent(out) :: P(1:ldp,1:n)

  call unif_dense(n, Q, ldq, P, ldp, qv, ufact)
end subroutine f90_unif_dense

subroutine f90_unif_csr(n, spQ, rowptr, colind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: rowptr(1:*), colind(1:*)

  call unif_csr(n, spQ, rowptr, colind, nnz, spP, qv, ufact)
end subroutine f90_unif_csr

subroutine f90_unif_csc(n, spQ, colptr, rowind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: colptr(1:*), rowind(1:*)

  call unif_csc(n, spQ, colptr, rowind, nnz, spP, qv, ufact)
end subroutine f90_unif_csc

subroutine f90_unif_coo(n, spQ, rowind, colind, nnz, spP, qv, ufact)
  use unif_matrix
  integer, intent(in) :: n, nnz
  double precision, intent(in) :: ufact
  double precision, intent(out) :: qv
  double precision, intent(in) :: spQ(1:nnz)
  double precision, intent(out) :: spP(1:nnz)
  integer, intent(in) :: rowind(1:*), colind(1:*)

  call unif_coo(n, spQ, rowind, colind, nnz, spP, qv, ufact)
end subroutine f90_unif_coo

!!!

subroutine f90_arnoldi_dense(trans, n, A, lda, x, incx, m, H, ldh, V, ldv, &
  beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:lda,1:n), x(1:incx,1:n), tol
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_dense(trans, n, A, lda, x, incx, m, H, ldh, V, ldv, &
    beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_dense

subroutine f90_arnoldi_csr(trans, n, A, rowptr, colind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_csr(trans, n, A, rowptr, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_csr

subroutine f90_arnoldi_csc(trans, n, A, colptr, rowind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_csc(trans, n, A, colptr, rowind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_csc

subroutine f90_arnoldi_coo(trans, n, A, rowind, colind, nnz, &
  x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
  use arnoldi
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, incx, m, ldh, ldv, ite
  double precision, intent(in) :: A(1:nnz), x(1:incx,1:n), tol
  integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
  double precision, intent(out) :: H(1:ldh,1:m), V(1:ldv,1:m)
  double precision, intent(out) :: beta, rnorm
  integer, intent(out) :: info

  call arnoldi_coo(trans, n, A, rowind, colind, nnz, &
    x, incx, m, H, ldh, V, ldv, beta, rnorm, tol, ite, info)
end subroutine f90_arnoldi_coo

!!!!!!!!!!!

subroutine f90_mpow_dense(trans, n, MA, lda, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, lda, lde, m
  double precision, intent(in) :: MA(1:lda,1:n)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_dense(trans, n, MA, lda, ME, lde, m, info)
end subroutine f90_mpow_dense

subroutine f90_mpow_csr(trans, n, spMA, rowptr, colind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: rowptr(1:*), colind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_csr(trans, n, spMA, rowptr, colind, nnz, ME, lde, m, info)
end subroutine f90_mpow_csr

subroutine f90_mpow_csc(trans, n, spMA, colptr, rowind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: colptr(1:*), rowind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_csc(trans, n, spMA, colptr, rowind, nnz, ME, lde, m, info)
end subroutine f90_mpow_csc

subroutine f90_mpow_coo(trans, n, spMA, rowind, colind, nnz, ME, lde, m, info)
  use mpow
  character, intent(in) :: trans
  integer, intent(in) :: n, nnz, lde, m
  integer, intent(in) :: rowind(1:*), colind(1:*)
  double precision, intent(in) :: spMA(1:*)
  double precision, intent(out) :: ME(1:lde,1:n)
  integer, intent(out) :: info

  call mpow_coo(trans, n, spMA, rowind, colind, nnz, ME, lde, m, info)
end subroutine f90_mpow_coo
