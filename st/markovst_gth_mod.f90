
module markovst_gth
  implicit none
  private
  public markovst_gth_dense
  public markovst_gth_csr
  public markovst_gth_csc
  public markovst_gth_coo

contains

  ! Description: stationary vector via GTH algorithm
  !        x * Q = 0, x*1 = 1
  !
  ! Q: CTMC kernel

  subroutine markovst_gth_base(n, A, lda, x, incx)
    integer, intent(in) :: n, lda, incx
    double precision, intent(inout) :: A(1:lda,1:n)
    double precision, intent(out) :: x(1:incx,1:n)

    integer :: i, j, l
    double precision :: tmp

    do l = n, 2, -1
       tmp = sum(A(l,1:(l-1)))
       do j = 1, l-1
          do i = 1, l-1
             if (i /= j) then
                A(i,j) = A(i,j) + A(l,j) * A(i,l) / tmp
             end if
          end do
       end do
       do i = l-1, 1, -1
          A(i,l) = A(i,l) / tmp
       end do
       do i = 1, l-1
          A(l,i) = 0.0d0
       end do
       A(l,l) = -1.0d0
    end do

    x(1,1) = 1.0d0
    do l = 2, n
       x(1,l) = sum(x(1,1:(l-1)) * A(1:(l-1),l))
    end do
    x(1,1:n) = x(1,1:n) / sum(x(1,1:n))
  end subroutine markovst_gth_base

  subroutine markovst_gth_dense(n, Q, ldq, x, incx)
    integer, intent(in) :: n, ldq, incx
    double precision, intent(in) :: Q(1:ldq,1:n)
    double precision, intent(out) :: x(1:incx,1:n)
    double precision :: A(1:n,1:n)

    A(1:n,1:n) = Q(1:n,1:n)
    call markovst_gth_base(n, A, n, x, incx)
  end subroutine markovst_gth_dense

  subroutine markovst_gth_csr(n, spQ, rowptr, colind, nnz, x, incx)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: rowptr(1:n+1), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)

    double precision :: A(1:n,1:n)

    call csr_to_dense(n, n, spQ, rowptr, colind, nnz, A, n)
    call markovst_gth_base(n, A, n, x, incx)
  end subroutine markovst_gth_csr

  subroutine markovst_gth_csc(n, spQ, colptr, rowind, nnz, x, incx)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: colptr(1:n+1), rowind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)

    double precision :: A(1:n,1:n)

    call csc_to_dense(n, n, spQ, colptr, rowind, nnz, A, n)
    call markovst_gth_base(n, A, n, x, incx)
  end subroutine markovst_gth_csc

  subroutine markovst_gth_coo(n, spQ, rowind, colind, nnz, x, incx)
    use sparse
    integer, intent(in) :: n, nnz, incx
    integer, intent(in) :: rowind(1:nnz), colind(1:nnz)
    double precision, intent(in) :: spQ(1:nnz)
    double precision, intent(out) :: x(1:incx,1:n)

    double precision :: A(1:n,1:n)

    call coo_to_dense(n, n, spQ, rowind, colind, nnz, A, n)
    call markovst_gth_base(n, A, n, x, incx)
  end subroutine markovst_gth_coo

end module markovst_gth

