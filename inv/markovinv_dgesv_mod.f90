
module markovinv_dgesv
  implicit none
  private
  public markovinv_dgesv_dense
  public markovinv_dgesv_csr
  public markovinv_dgesv_csc
  public markovinv_dgesv_coo

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

  subroutine markovinv_dgesv_base(trans, n, nrhs, A, lda, x, ldx, y, ldy, info)
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, lda, ldx, ldy
    double precision, intent(inout) :: A(1:lda,1:n)
    double precision, intent(in) :: x(1:ldx,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs)
    integer, intent(out) :: info

    integer :: ipiv(1:n)

    select case (trans)
      case ('N', 'n')
        A(1:n,1:n) = -A(1:n,1:n)
      case ('T', 't')
        A(1:n,1:n) = -transpose(A(1:n,1:n))
    end select

    y(1:n,1:nrhs) = x(1:n,1:nrhs)
    call dgesv(n, nrhs, A, lda, ipiv, y, ldy, info)
  end subroutine markovinv_dgesv_base

  subroutine markovinv_dgesv_dense(trans, n, nrhs, Q, ldq, x, ldx, y, ldy, info)
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, ldq, ldx, ldy
    double precision, intent(in) :: Q(1:ldq,1:n)
    double precision, intent(in) :: x(1:ldx,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    A(1:n,1:n) = Q(1:n,1:n)
    call markovinv_dgesv_base(trans, n, nrhs, A, n, x, ldx, y, ldy, info)
  end subroutine markovinv_dgesv_dense

  subroutine markovinv_dgesv_csr(trans, n, nrhs, spQ, rowptr, colind, nnz, x, ldx, y, ldy, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, nnz, ldx, ldy
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(in) :: x(1:ldx,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csr_to_dense(n, n, spQ, rowptr, colind, nnz, A, n)
    call markovinv_dgesv_base(trans, n, nrhs, A, n, x, ldx, y, ldy, info)
  end subroutine markovinv_dgesv_csr

  subroutine markovinv_dgesv_csc(trans, n, nrhs, spQ, colptr, rowind, nnz, x, ldx, y, ldy, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, nnz, ldx, ldy
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csc_to_dense(n, n, spQ, colptr, rowind, nnz, A, n)
    call markovinv_dgesv_base(trans, n, nrhs, A, n, x, ldx, y, ldy, info)
  end subroutine markovinv_dgesv_csc

  subroutine markovinv_dgesv_coo(trans, n, nrhs, spQ, rowind, colind, nnz, x, ldx, y, ldy, info)
    use sparse
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, nnz, ldx, ldy
    integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz), x(1:ldx,1:nrhs)
    double precision, intent(out) :: y(1:ldy,1:nrhs)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call coo_to_dense(n, n, spQ, rowind, colind, nnz, A, n)
    call markovinv_dgesv_base(trans, n, nrhs, A, n, x, ldx, y, ldy, info)
  end subroutine markovinv_dgesv_coo

end module markovinv_dgesv

