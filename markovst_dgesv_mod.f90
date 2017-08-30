
module markovst_dgesv
  implicit none
  private
  public markovst_dgesv_dense
  public markovst_dgesv_csr
  public markovst_dgesv_csc
  public markovst_dgesv_coo

contains

  ! Description: stationary vector via GTH algorithm
  !        x * Q = 0, x*1 = 1
  !
  ! Q: CTMC kernel

  subroutine markovst_dgesv_base(n, A, lda, x, incx, info)
    integer, parameter :: index = 1
    double precision, parameter :: value = 1.0d0
    integer, intent(in) :: n, lda, incx
    double precision, intent(inout) :: A(1:lda,1:n)
    double precision, intent(out) :: x(1:incx,1:n)
    integer, intent(out) :: info

    integer :: ipiv(1:n)
    double precision :: b(1:n)

    A(1:n,index) = value
    A(1:n,1:n) = transpose(A(1:n,1:n))
    b(1:n) = 0.0d0
    b(index) = value
    call dgesv(n, 1, A, lda, ipiv, b, n, info)
    call dcopy(n, b, 1, x, incx)
  end subroutine markovst_dgesv_base

  subroutine markovst_dgesv_dense(n, Q, ldq, x, incx, info)
    integer, intent(in) :: n, ldq, incx
    double precision, intent(in) :: Q(1:ldq,1:n)
    double precision, intent(out) :: x(1:incx,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    A(1:n,1:n) = Q(1:n,1:n)
    call markovst_dgesv_base(n, A, n, x, incx, info)
  end subroutine markovst_dgesv_dense

  subroutine markovst_dgesv_csr(n, spQ, rowptr, colind, nnz, x, incx, info)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csr_to_dense(n, n, spQ, rowptr, colind, nnz, A, n)
    call markovst_dgesv_base(n, A, n, x, incx, info)
  end subroutine markovst_dgesv_csr

  subroutine markovst_dgesv_csc(n, spQ, colptr, rowind, nnz, x, incx, info)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csc_to_dense(n, n, spQ, colptr, rowind, nnz, A, n)
    call markovst_dgesv_base(n, A, n, x, incx, info)
  end subroutine markovst_dgesv_csc

  subroutine markovst_dgesv_coo(n, spQ, rowind, colind, nnz, x, incx, info)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call coo_to_dense(n, n, spQ, rowind, colind, nnz, A, n)
    call markovst_dgesv_base(n, A, n, x, incx, info)
  end subroutine markovst_dgesv_coo

end module markovst_dgesv

