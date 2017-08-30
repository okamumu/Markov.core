!
! wrappers for markovinv_gs_mod.f90
!

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
