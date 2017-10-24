
module markovstsen_dgesv
  implicit none
  private
  public markovstsen_dgesv_dense
  public markovstsen_dgesv_csr
  public markovstsen_dgesv_csc
  public markovstsen_dgesv_coo

contains

  ! Description: sensitivity vector for stationary vector
  !
  !        s * Q + pis * dQ = 0,  s*1 = 0
  !
  ! Q: CTMC kernel
  ! dQ: the first derivative of CTMC kernel with respect to a parameter
  ! pis: stationary vector such as pis * Q = 0, pis*1 = 1

  subroutine markovstsen_dgesv_base(n, A, lda, s, incs, info)
    integer, parameter :: index = 1
    double precision, parameter :: value = 1.0d0
    integer, intent(in) :: n, lda, incs
    double precision, intent(inout) :: A(1:lda,1:n)
    double precision, intent(out) :: s(1:incs,1:n)
    integer, intent(out) :: info

    integer :: ipiv(1:n)
    double precision :: b(1:n)

    A(1:n,index) = value
    A(1:n,1:n) = transpose(A(1:n,1:n))
    call dcopy(n, s, incs, b, 1)
    b(index) = 0.0d0
    call dgesv(n, 1, A, lda, ipiv, b, n, info)
    call dcopy(n, b, 1, s, incs)
  end subroutine markovstsen_dgesv_base

  subroutine markovstsen_dgesv_dense(n, Q, ldq, dQ, lddq, pis, incp, s, incs, info)
    integer, intent(in) :: n, ldq, lddq, incp, incs
    double precision, intent(in) :: Q(1:ldq,1:n), dQ(1:lddq,1:n)
    double precision, intent(in) :: pis(1:incp,1:n)
    double precision, intent(out) :: s(1:incs,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    A(1:n,1:n) = Q(1:n,1:n)
    call dgemv('T', n, n, -1.0d0, dQ, lddq, pis, incp, 0.0d0, s, incs)
    call markovstsen_dgesv_base(n, A, n, s, incs, info)
  end subroutine markovstsen_dgesv_dense

  subroutine markovstsen_dgesv_csr(n, spQ, rowptr0, colind0, nnz0, &
    dQ, rowptr1, colind1, nnz1, pis, incp, s, incs, info)
    use sparse
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, incp, incs
    integer, intent(in) :: rowptr0(1:n+1), colind0(1:nnz0)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: spQ(1:nnz0)
    double precision, intent(in) :: dQ(1:nnz1)
    double precision, intent(in) :: pis(1:incp,1:n)
    double precision, intent(out) :: s(1:incs,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csr_to_dense(n, n, spQ, rowptr0, colind0, nnz0, A, n)
    call dcsrmv('T', n, n, -1.0d0, dQ, rowptr1, colind1, nnz1, pis, incp, 0.0d0, s, incs)
    call markovstsen_dgesv_base(n, A, n, s, incs, info)
  end subroutine markovstsen_dgesv_csr

  subroutine markovstsen_dgesv_csc(n, spQ, colptr0, rowind0, nnz0, &
    dQ, colptr1, rowind1, nnz1, pis, incp, s, incs, info)
    use sparse
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, incp, incs
    integer, intent(in) :: colptr0(1:n+1), rowind0(1:nnz0)
    integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
    double precision, intent(in) :: spQ(1:nnz0)
    double precision, intent(in) :: dQ(1:nnz1)
    double precision, intent(in) :: pis(1:incp,1:n)
    double precision, intent(out) :: s(1:incs,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call csc_to_dense(n, n, spQ, colptr0, rowind0, nnz0, A, n)
    call dcscmv('T', n, n, -1.0d0, dQ, colptr1, rowind1, nnz1, pis, incp, 0.0d0, s, incs)
    call markovstsen_dgesv_base(n, A, n, s, incs, info)
  end subroutine markovstsen_dgesv_csc

  subroutine markovstsen_dgesv_coo(n, spQ, rowind0, colind0, nnz0, &
    dQ, rowind1, colind1, nnz1, pis, incp, s, incs, info)
    use sparse
    use spblas
    integer, intent(in) :: n, nnz0, nnz1, incp, incs
    integer, intent(in) :: rowind0(1:nnz0), colind0(1:nnz0)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: spQ(1:nnz0)
    double precision, intent(in) :: dQ(1:nnz1)
    double precision, intent(in) :: pis(1:incp,1:n)
    double precision, intent(out) :: s(1:incs,1:n)
    integer, intent(out) :: info

    double precision :: A(1:n,1:n)

    call coo_to_dense(n, n, spQ, rowind0, colind0, nnz0, A, n)
    call dcoomv('T', n, n, -1.0d0, dQ, rowind1, colind1, nnz1, pis, incp, 0.0d0, s, incs)
    call markovstsen_dgesv_base(n, A, n, s, incs, info)
  end subroutine markovstsen_dgesv_coo

end module markovstsen_dgesv

