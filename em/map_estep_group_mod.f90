
module map_estep_group
  implicit none
  private
  public map_estep_group_unif_dense
  public map_estep_group_unif_csr
  public map_estep_group_unif_csc
  public map_estep_group_unif_coo

contains

! Description: estep for PH with weighted time and group/truncated data

! alpha      (in): initial vector
! baralpha   (in): baralpha = alpha (-Q)^-1
! xi         (in): exit vector
! one        (in): one vector
! Q          (in): infinitesimal generator
! P          (in): uniformed generator
! qv         (in): uniformization constant
! tdat       (in): interarrival time
! wdat       (in): weights for interarrivals
! gdat       (in): # of arrivals (-1 means NA)
! gdatlast   (in): # of arrivals in [lasttime, infinity] (-1 means NA)
! idat       (in): indicator whether an arrival occurs at the last instant
! etotal    (out): expected # of arrivals
! eb        (out): expected # of starts
! ey        (out): expected # of exits
! ez        (out): expected sojourn time
! en        (out): expected # of phase transitions
! poi_eps    (in): eps for poisson prob (optional)
! atol       (in): tolerance error in uniformization (optional)

! return value -> llf (log-likelihood)

  subroutine map_estep_group_unif_dense(n, alpha, ev, &
    D0, P0, D1, P1, qv, poieps, &
    m, tdat, gdat, idat, &
    eb, ez, en0, en1, llf)
    use poisson
    use map_mexp_unif
    use map_mexpconv_unif
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m
    double precision, intent(in) :: alpha(1:n), ev(1:n)
    double precision, intent(in) :: D0(1:n,1:n), D1(1:n,1:n)
    double precision, intent(in) :: P0(1:n,1:n), P1(1:n,1:n)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), idat(1:m)
    double precision, intent(out) :: eb(1:n), ez(1:n)
    double precision, intent(out) :: en0(1:n,1:n), en1(1:n,1:n)
    double precision, intent(out) :: llf

    integer :: k, right, nmax
    double precision :: weight, scale, tmax
    double precision :: ddot

    double precision :: tmpv(1:n), tmpz(1:n)
    double precision :: hen0(1:n,1:n), hen1(1:n,1:n)
    double precision :: vf(1:n,0:m+1), vb(1:n,0:m+1)

    double precision, allocatable :: poi(:)
    double precision, allocatable :: xi(:,:), vc(:,:,:)

    eb = 0.0d0
    ez = 0.0d0
    en0 = 0.0d0
    en1 = 0.0d0
    llf = 0.0d0

    tmax = maxval(tdat(1:m))
    nmax = maxval(gdat(1:m))
    right = max(poisson_rightbound(qv*tmax, poieps), nmax) + 1
    allocate(poi(0:right))
    allocate(xi(1:n,0:nmax))
    allocate(vc(1:n,0:nmax,left:right))

    ! backward
    call dcopy(n, ev, 1, vb(1,m+1), 1)
    do k = m, 1, -1
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dgemv('N', n, n, 1.0d0, D1, n, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      call map_mexp_unif_dense('N', n, P0, n, P1, n, qv, &
        0, right, poi, weight, gdat(k), tmpv, 1, vb(1,k), 1, atol, xi)

      scale = sum(vb(1:n,k))
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      llf = llf + log(scale)
    end do
    call dcopy(n, vb(1,1), 1, eb, 1)

    ! forward
    call dcopy(n, alpha, 1, vf(1,0), 1)
    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call map_mexp_unif_dense('T', n, P0, n, P1, n, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, atol, xi)

      if (idat(k) == 1) then
        call dgemv('T', n, n, 1.0d0, D1, n, tmpv, 1, 0.0d0, vb(1,k), 1)
      else
        call dcopy(n, tmpv, 1, vf(1,k), 1)
      end if

      scale = sum(vf(1:n,k))
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
    end do

    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dgemv('N', n, n, 1.0d0, D1, n, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      hen0 = 0.0d0
      hen1 = 0.0d0
      call map_mexpconv_unif_dense('T', 'N', n, P0, n, P1, n, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, tmpz, 1, &
        hen0, n, hen1, n, atol, xi, vc)
      scale = ddot(n, tmpz, 1, tmpv, 1)
      call daxpy(n*n, 1.0d0/scale, hen0, 1, en0, 1)
      call daxpy(n*n, 1.0d0/scale, hen1, 1, en1, 1)

      if (idat(k) == 1) then
        call dger(n, n, 1.0/scale, tmpz, 1, vb(1,k+1), 1, en1, n)
      end if
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    scale = sum(eb(1:n))
    call dscal(n, 1.0d0/scale, eb, 1)
    llf = llf + log(scale)
    call dcopy(n, en0(1,1), n+1, ez, 1)
    en0(1:n,1:n) = D0(1:n,1:n) * en0(1:n,1:n)
    en1(1:n,1:n) = D1(1:n,1:n) * en1(1:n,1:n)

    deallocate(poi, xi, vc)
  end subroutine map_estep_group_unif_dense

  subroutine map_estep_group_unif_csr(n, alpha, ev, &
    D0, P0, rowptr0, colind0, nnz0, &
    D1, P1, rowptr1, colind1, nnz1, &
    qv, poieps, &
    m, tdat, gdat, idat, &
    eb, ez, en0, en1, llf)
    use sparse
    use spblas
    use poisson
    use map_mexp_unif
    use map_mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m, nnz0, nnz1
    double precision, intent(in) :: alpha(1:n), ev(1:n)
    double precision, intent(in) :: D0(1:nnz0), D1(1:nnz1)
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowptr0(base:base+n), colind0(base:base+nnz0-1)
    integer, intent(in) :: rowptr1(1:n+1), colind1(1:nnz1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), idat(1:m)
    double precision, intent(out) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en0(base:base+nnz0-1), en1(1:nnz1)
    double precision, intent(out) :: llf

    integer :: i, j, z, k, right, nmax
    double precision :: weight, scale, tmax
    double precision :: ddot

    double precision :: tmpv(1:n), tmpz(1:n)
    double precision :: hen0(1:nnz0), hen1(1:nnz1)
    double precision :: vf(1:n,0:m+1), vb(1:n,0:m+1)

    double precision, allocatable :: poi(:)
    double precision, allocatable :: xi(:,:), vc(:,:,:)

    eb = 0.0d0
    ez = 0.0d0
    en0 = 0.0d0
    en1 = 0.0d0
    llf = 0.0d0

    tmax = maxval(tdat(1:m))
    nmax = maxval(gdat(1:m))
    right = max(poisson_rightbound(qv*tmax, poieps), nmax) + 1
    allocate(poi(0:right))
    allocate(xi(1:n,0:nmax))
    allocate(vc(1:n,0:nmax,left:right))

    ! backward
    call dcopy(n, ev, 1, vb(1,m+1), 1)
    do k = m, 1, -1
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcsrmv('N', n, n, 1.0d0, D1, rowptr1, colind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      call map_mexp_unif_csr('N', n, P0, rowptr0, colind0, nnz0, &
        P1, rowptr1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), tmpv, 1, vb(1,k), 1, atol, xi)

      scale = sum(vb(1:n,k))
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      llf = llf + log(scale)
    end do
    call dcopy(n, vb(1,1), 1, eb, 1)

    ! forward
    call dcopy(n, alpha, 1, vf(1,0), 1)
    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call map_mexp_unif_csr('T', n, P0, rowptr0, colind0, nnz0, &
        P1, rowptr1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, atol, xi)

      if (idat(k) == 1) then
        call dcsrmv('T', n, n, 1.0d0, D1, rowptr1, colind1, nnz1, tmpv, 1, 0.0d0, vb(1,k), 1)
      else
        call dcopy(n, tmpv, 1, vf(1,k), 1)
      end if

      scale = sum(vf(1:n,k))
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
    end do

    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcsrmv('N', n, n, 1.0d0, D1, rowptr1, colind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      hen0 = 0.0d0
      hen1 = 0.0d0
      call map_mexpconv_unif_csr('T', 'N', n, P0, rowptr0, colind0, nnz0, &
        P1, rowptr1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, tmpz, 1, &
        hen0, hen1, atol, xi, vc)
      scale = ddot(n, tmpz, 1, tmpv, 1)
      call daxpy(nnz0, 1.0d0/scale, hen0, 1, en0, 1)
      call daxpy(nnz1, 1.0d0/scale, hen1, 1, en1, 1)

      if (idat(k) == 1) then
        call dcsrr(n, n, 1.0/scale, tmpz, 1, vb(1,k+1), 1, en1, rowptr1, colind1, nnz1)
      end if
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    scale = sum(eb(1:n))
    call dscal(n, 1.0d0/scale, eb, 1)
    llf = llf + log(scale)
    do i = base, base+n-1
      do z = rowptr0(i), rowptr0(i+1)-1
        j = colind0(z)
        if (i == j) then
          ez(i) = en0(z)
          exit
        end if
      end do
    end do
    en0(:) = D0(:) * en0(:)
    en1(:) = D1(:) * en1(:)

    deallocate(poi, xi, vc)
  end subroutine map_estep_group_unif_csr

  subroutine map_estep_group_unif_csc(n, alpha, ev, &
    D0, P0, colptr0, rowind0, nnz0, &
    D1, P1, colptr1, rowind1, nnz1, &
    qv, poieps, &
    m, tdat, gdat, idat, &
    eb, ez, en0, en1, llf)
    use sparse
    use spblas
    use poisson
    use map_mexp_unif
    use map_mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m, nnz0, nnz1
    double precision, intent(in) :: alpha(1:n), ev(1:n)
    double precision, intent(in) :: D0(1:nnz0), D1(1:nnz1)
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: colptr0(base:base+n), rowind0(base:base+nnz0-1)
    integer, intent(in) :: colptr1(1:n+1), rowind1(1:nnz1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), idat(1:m)
    double precision, intent(out) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en0(base:base+nnz0-1), en1(1:nnz1)
    double precision, intent(out) :: llf

    integer :: i, j, z, k, right, nmax
    double precision :: weight, scale, tmax
    double precision :: ddot

    double precision :: tmpv(1:n), tmpz(1:n)
    double precision :: hen0(1:nnz0), hen1(1:nnz1)
    double precision :: vf(1:n,0:m+1), vb(1:n,0:m+1)

    double precision, allocatable :: poi(:)
    double precision, allocatable :: xi(:,:), vc(:,:,:)

    eb = 0.0d0
    ez = 0.0d0
    en0 = 0.0d0
    en1 = 0.0d0
    llf = 0.0d0

    tmax = maxval(tdat(1:m))
    nmax = maxval(gdat(1:m))
    right = max(poisson_rightbound(qv*tmax, poieps), nmax) + 1
    allocate(poi(0:right))
    allocate(xi(1:n,0:nmax))
    allocate(vc(1:n,0:nmax,left:right))

    ! backward
    call dcopy(n, ev, 1, vb(1,m+1), 1)
    do k = m, 1, -1
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcscmv('N', n, n, 1.0d0, D1, colptr1, rowind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      call map_mexp_unif_csc('N', n, P0, colptr0, rowind0, nnz0, &
        P1, colptr1, rowind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), tmpv, 1, vb(1,k), 1, atol, xi)

      scale = sum(vb(1:n,k))
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      llf = llf + log(scale)
    end do
    call dcopy(n, vb(1,1), 1, eb, 1)

    ! forward
    call dcopy(n, alpha, 1, vf(1,0), 1)
    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call map_mexp_unif_csc('T', n, P0, colptr0, rowind0, nnz0, &
        P1, colptr1, rowind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, atol, xi)

      if (idat(k) == 1) then
        call dcscmv('T', n, n, 1.0d0, D1, colptr1, rowind1, nnz1, tmpv, 1, 0.0d0, vb(1,k), 1)
      else
        call dcopy(n, tmpv, 1, vf(1,k), 1)
      end if

      scale = sum(vf(1:n,k))
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
    end do

    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcscmv('N', n, n, 1.0d0, D1, colptr1, rowind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      hen0 = 0.0d0
      hen1 = 0.0d0
      call map_mexpconv_unif_csc('T', 'N', n, P0, colptr0, rowind0, nnz0, &
        P1, colptr1, rowind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, tmpz, 1, &
        hen0, hen1, atol, xi, vc)
      scale = ddot(n, tmpz, 1, tmpv, 1)
      call daxpy(nnz0, 1.0d0/scale, hen0, 1, en0, 1)
      call daxpy(nnz1, 1.0d0/scale, hen1, 1, en1, 1)

      if (idat(k) == 1) then
        call dcscr(n, n, 1.0/scale, tmpz, 1, vb(1,k+1), 1, en1, colptr1, rowind1, nnz1)
      end if
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    scale = sum(eb(1:n))
    call dscal(n, 1.0d0/scale, eb, 1)
    llf = llf + log(scale)
    do j = base, base+n-1
      do z = colptr0(j), colptr0(j+1)-1
        i = rowind0(z)
        if (i == j) then
          ez(i) = en0(z)
          exit
        end if
      end do
    end do
    en0(:) = D0(:) * en0(:)
    en1(:) = D1(:) * en1(:)

    deallocate(poi, xi, vc)
  end subroutine map_estep_group_unif_csc

  subroutine map_estep_group_unif_coo(n, alpha, ev, &
    D0, P0, rowind0, colind0, nnz0, &
    D1, P1, rowind1, colind1, nnz1, &
    qv, poieps, &
    m, tdat, gdat, idat, &
    eb, ez, en0, en1, llf)
    use sparse
    use spblas
    use poisson
    use map_mexp_unif
    use map_mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m, nnz0, nnz1
    double precision, intent(in) :: alpha(1:n), ev(1:n)
    double precision, intent(in) :: D0(1:nnz0), D1(1:nnz1)
    double precision, intent(in) :: P0(1:nnz0), P1(1:nnz1)
    integer, intent(in) :: rowind0(base:base+nnz0-1), colind0(base:base+nnz0-1)
    integer, intent(in) :: rowind1(1:nnz1), colind1(1:nnz1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), idat(1:m)
    double precision, intent(out) :: eb(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en0(base:base+nnz0-1), en1(1:nnz1)
    double precision, intent(out) :: llf

    integer :: i, j, z, k, right, nmax
    double precision :: weight, scale, tmax
    double precision :: ddot

    double precision :: tmpv(1:n), tmpz(1:n)
    double precision :: hen0(1:nnz0), hen1(1:nnz1)
    double precision :: vf(1:n,0:m+1), vb(1:n,0:m+1)

    double precision, allocatable :: poi(:)
    double precision, allocatable :: xi(:,:), vc(:,:,:)

    eb = 0.0d0
    ez = 0.0d0
    en0 = 0.0d0
    en1 = 0.0d0
    llf = 0.0d0

    tmax = maxval(tdat(1:m))
    nmax = maxval(gdat(1:m))
    right = max(poisson_rightbound(qv*tmax, poieps), nmax) + 1
    allocate(poi(0:right))
    allocate(xi(1:n,0:nmax))
    allocate(vc(1:n,0:nmax,left:right))

    ! backward
    call dcopy(n, ev, 1, vb(1,m+1), 1)
    do k = m, 1, -1
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcoomv('N', n, n, 1.0d0, D1, rowind1, colind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      call map_mexp_unif_coo('N', n, P0, rowind0, colind0, nnz0, &
        P1, rowind1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), tmpv, 1, vb(1,k), 1, atol, xi)

      scale = sum(vb(1:n,k))
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      llf = llf + log(scale)
    end do
    call dcopy(n, vb(1,1), 1, eb, 1)

    ! forward
    call dcopy(n, alpha, 1, vf(1,0), 1)
    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call map_mexp_unif_coo('T', n, P0, rowind0, colind0, nnz0, &
        P1, rowind1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, atol, xi)

      if (idat(k) == 1) then
        call dcoomv('T', n, n, 1.0d0, D1, rowind1, colind1, nnz1, tmpv, 1, 0.0d0, vb(1,k), 1)
      else
        call dcopy(n, tmpv, 1, vf(1,k), 1)
      end if

      scale = sum(vf(1:n,k))
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
    end do

    do k = 1, m
      right = max(poisson_rightbound(qv*tdat(k), poieps), gdat(k)) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      if (idat(k) == 1) then
        call dcoomv('N', n, n, 1.0d0, D1, rowind1, colind1, nnz1, vb(1,k+1), 1, 0.0d0, tmpv, 1)
      else
        call dcopy(n, vb(1,k+1), 1, tmpv, 1)
      end if

      hen0 = 0.0d0
      hen1 = 0.0d0
      call map_mexpconv_unif_coo('T', 'N', n, P0, rowind0, colind0, nnz0, &
        P1, rowind1, colind1, nnz1, qv, &
        0, right, poi, weight, gdat(k), vf(1,k-1), 1, tmpv, 1, tmpz, 1, &
        hen0, hen1, atol, xi, vc)
      scale = ddot(n, tmpz, 1, tmpv, 1)
      call daxpy(nnz0, 1.0d0/scale, hen0, 1, en0, 1)
      call daxpy(nnz1, 1.0d0/scale, hen1, 1, en1, 1)

      if (idat(k) == 1) then
        call dcoor(n, n, 1.0/scale, tmpz, 1, vb(1,k+1), 1, en1, rowind1, colind1, nnz1)
      end if
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    scale = sum(eb(1:n))
    call dscal(n, 1.0d0/scale, eb, 1)
    llf = llf + log(scale)
    do z = base, base+nnz0-1
      i = rowind0(z)
      j = colind0(z)
      if (i == j) then
        ez(i) = en0(z)
      end if
    end do
    en0(:) = D0(:) * en0(:)
    en1(:) = D1(:) * en1(:)

    deallocate(poi, xi, vc)
  end subroutine map_estep_group_unif_coo

end module map_estep_group

