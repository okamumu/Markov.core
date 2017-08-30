!
! gph estep
!

module gph_estep_grouppoi
  implicit none

contains

! Description: estep for PH with weighted time and group/truncated data

! alpha      (in): initial vector
! baralpha   (in): baralpha = alpha (-Q)^-1
! xi         (in): exit vector
! ev         (in): one vector
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

  subroutine gph_estep_grouppoi_unif_dense(n, alpha, baralpha, tau, ev, omega, &
    T, P, qv, poieps, m, tdat, gdat, gdatlast, idat, &
    etotal, eb, ey, ez, en, llf)
    use gamma
    use poisson
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m
    double precision, intent(in) :: omega
    double precision, intent(in) :: alpha(1:n), baralpha(1:n), tau(1:n), ev(1:n)
    double precision, intent(in) :: T(1:n,1:n), P(1:n,1:n), qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(out) :: en(1:n,1:n)
    double precision, intent(out) :: llf

    integer :: k, right
    double precision :: weight, tmax
    double precision :: tmp, nn, uu
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: tildevf(1:n), tildevb(1:n)
    double precision :: barvf(1:n,0:m), barvb(1:n,0:m)
    double precision :: vf(1:n), vb(1:n,0:m)
    double precision :: vc(1:n,0:m)
    double precision :: wg(1:m+1), wp(1:m)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, baralpha, 1, barvf(1,0), 1)
    call dcopy(n, ev, 1, barvb(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_dense_vec('T', n, P, n, qv, left, right, poi, weight, &
        barvf(1,k-1), 1, barvf(1,k), 1, atol)
      call mexp_unif_dense_vec('N', n, P, n, qv, left, right, poi, weight, &
        barvb(1,k-1), 1, barvb(1,k), 1, atol)

      call dgemv('N', n, n, -1.0d0, T, n, barvb(1,k), 1, 0.0d0, vb(1,k), 1)

      tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
      tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)

      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmp = ddot(n, alpha, 1, tildevb, 1)
        llf = llf + gdat(k) * log(tmp) - logfact(gdat(k))
        nn = nn + gdat(k)
        uu = uu + tmp
        wg(k) = gdat(k) / tmp
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if

      if (idat(k) == 1) then
        call dgemv('T', n, n, -1.0d0, T, n, barvf(1,k), 1, 0.0d0, vf, 1)
        tmp = ddot(n, alpha, 1, vb(1,k), 1)
        llf = llf + log(tmp)
        nn = nn + 1.0d0
        wp(k) = 1.0d0 / tmp
        call daxpy(n, wp(k), vb(1,k), 1, eb, 1)
        call daxpy(n, wp(k), vf, 1, ey, 1)
      end if
    end do

    ! for the interval [t_m, infinity)
    if (gdatlast >= 0) then
      tmp = ddot(n, alpha, 1, barvb(1,m), 1)
      llf = llf + gdatlast * log(tmp) - logfact(gdatlast)
      nn = nn + gdatlast
      uu = uu + tmp
      wg(m+1) = gdatlast / tmp
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if

    ! compupte weights for unobserved periods
    do k = 1, m
      if (gdat(k) == -1) then
        tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
        tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)
        wg(k) = omega
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if
    end do
    if (gdatlast == -1) then
      wg(m+1) = omega
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if
    llf = llf + nn * log(omega) - omega * uu

    ! compute vectors for convolution
    vc(1:n,m) = 0.0d0
    call daxpy(n, wg(m+1)-wg(m), baralpha, 1, vc(1,m), 1)
    if (idat(m) == 1) then
      call daxpy(n, wp(m), alpha, 1, vc(1,m), 1)
    end if
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_dense_vec('T', n, P, n, qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call daxpy(n, wg(k+1)-wg(k), baralpha, 1, vc(1,k), 1)

      if (idat(k) == 1) then
        call daxpy(n, wp(k), alpha, 1, vc(1,k), 1)
      end if
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call dger(n, n, wg(k+1)-wg(k), baralpha, 1, barvb(1,k), 1, en, n)
      call mexpconv_unif_dense_vec('T', 'N', n, P, n, qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, n, atol)
    end do
    call dger(n, n, wg(1), baralpha, 1, barvb(1,0), 1, en, n)

    etotal = nn + omega * (1.0d0 - uu)
    eb(1:n) = alpha(1:n) * eb(1:n)
    ey(1:n) = tau(1:n) * ey(1:n)
    call dcopy(n, en(1,1), n+1, ez, 1)
    en(1:n,1:n) = T(1:n,1:n) * en(1:n,1:n)

    deallocate(poi)
  end subroutine gph_estep_grouppoi_unif_dense

  subroutine gph_estep_grouppoi_unif_csr(n, alpha, baralpha, tau, ev, omega, &
    spT, spP, rowptr, colind, nnz, qv, poieps, &
    m, tdat, gdat, gdatlast, idat, &
    etotal, eb, ey, ez, en, llf)
    use sparse
    use spblas
    use gamma
    use poisson
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: omega
    double precision, intent(in) :: alpha(1:n), baralpha(1:n), tau(1:n), ev(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    integer, intent(in) :: rowptr(base:base+n), colind(base:base+nnz-1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, tmax
    double precision :: tmp, nn, uu
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: tildevf(1:n), tildevb(1:n)
    double precision :: barvf(1:n,0:m), barvb(1:n,0:m)
    double precision :: vf(1:n), vb(1:n,0:m)
    double precision :: vc(1:n,0:m)
    double precision :: wg(1:m+1), wp(1:m)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, baralpha, 1, barvf(1,0), 1)
    call dcopy(n, ev, 1, barvb(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_csr_vec('T', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        barvf(1,k-1), 1, barvf(1,k), 1, atol)
      call mexp_unif_csr_vec('N', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        barvb(1,k-1), 1, barvb(1,k), 1, atol)

      call dcsrmv('N', n, n, -1.0d0, spT, rowptr, colind, nnz, &
        barvb(1,k), 1, 0.0d0, vb(1,k), 1)

      tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
      tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)

      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmp = ddot(n, alpha, 1, tildevb, 1)
        llf = llf + gdat(k) * log(tmp) - logfact(gdat(k))
        nn = nn + gdat(k)
        uu = uu + tmp
        wg(k) = gdat(k) / tmp
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if

      if (idat(k) == 1) then
        call dcsrmv('T', n, n, -1.0d0, spT, rowptr, colind, nnz, &
          barvf(1,k), 1, 0.0d0, vf, 1)
        tmp = ddot(n, alpha, 1, vb(1,k), 1)
        llf = llf + log(tmp)
        nn = nn + 1.0d0
        wp(k) = 1.0d0 / tmp
        call daxpy(n, wp(k), vb(1,k), 1, eb, 1)
        call daxpy(n, wp(k), vf, 1, ey, 1)
      end if
    end do

    ! for the interval [t_m, infinity)
    if (gdatlast >= 0) then
      tmp = ddot(n, alpha, 1, barvb(1,m), 1)
      llf = llf + gdatlast * log(tmp) - logfact(gdatlast)
      nn = nn + gdatlast
      uu = uu + tmp
      wg(m+1) = gdatlast / tmp
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if

    ! compupte weights for unobserved periods
    do k = 1, m
      if (gdat(k) == -1) then
        tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
        tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)
        wg(k) = omega
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if
    end do
    if (gdatlast == -1) then
      wg(m+1) = omega
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if
    llf = llf + nn * log(omega) - omega * uu

    ! compute vectors for convolution
    vc(1:n,m) = 0.0d0
    call daxpy(n, wg(m+1)-wg(m), baralpha, 1, vc(1,m), 1)
    if (idat(m) == 1) then
      call daxpy(n, wp(m), alpha, 1, vc(1,m), 1)
    end if
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_csr_vec('T', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call daxpy(n, wg(k+1)-wg(k), baralpha, 1, vc(1,k), 1)

      if (idat(k) == 1) then
        call daxpy(n, wp(k), alpha, 1, vc(1,k), 1)
      end if
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call dcsrr(n, n, wg(k+1)-wg(k), baralpha, 1, barvb(1,k), 1, &
        en, rowptr, colind, nnz)
      call mexpconv_unif_csr_vec('T', 'N', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do
    call dcsrr(n, n, wg(1), baralpha, 1, barvb(1,0), 1, &
      en, rowptr, colind, nnz)

    etotal = nn + omega * (1.0d0 - uu)
    eb(1:n) = alpha(1:n) * eb(1:n)
    ey(1:n) = tau(1:n) * ey(1:n)
    do i = base, base+n-1
      do z = rowptr(i), rowptr(i+1)-1
        j = colind(z)
        if (i == j) then
          ez(i) = en(z)
          exit
        end if
      end do
    end do
    en(:) = spT(:) * en(:)

    deallocate(poi)
  end subroutine gph_estep_grouppoi_unif_csr

  subroutine gph_estep_grouppoi_unif_csc(n, alpha, baralpha, tau, ev, omega, &
    spT, spP, colptr, rowind, nnz, qv, poieps, &
    m, tdat, gdat, gdatlast, idat, &
    etotal, eb, ey, ez, en, llf)
    use sparse
    use spblas
    use gamma
    use poisson
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: omega
    double precision, intent(in) :: alpha(1:n), baralpha(1:n), tau(1:n), ev(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    integer, intent(in) :: colptr(base:base+n), rowind(base:base+nnz-1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, tmax
    double precision :: tmp, nn, uu
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: tildevf(1:n), tildevb(1:n)
    double precision :: barvf(1:n,0:m), barvb(1:n,0:m)
    double precision :: vf(1:n), vb(1:n,0:m)
    double precision :: vc(1:n,0:m)
    double precision :: wg(1:m+1), wp(1:m)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, baralpha, 1, barvf(1,0), 1)
    call dcopy(n, ev, 1, barvb(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_csc_vec('T', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        barvf(1,k-1), 1, barvf(1,k), 1, atol)
      call mexp_unif_csc_vec('N', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        barvb(1,k-1), 1, barvb(1,k), 1, atol)

      call dcscmv('N', n, n, -1.0d0, spT, colptr, rowind, nnz, &
        barvb(1,k), 1, 0.0d0, vb(1,k), 1)

      tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
      tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)

      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmp = ddot(n, alpha, 1, tildevb, 1)
        llf = llf + gdat(k) * log(tmp) - logfact(gdat(k))
        nn = nn + gdat(k)
        uu = uu + tmp
        wg(k) = gdat(k) / tmp
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if

      if (idat(k) == 1) then
        call dcscmv('T', n, n, -1.0d0, spT, colptr, rowind, nnz, &
          barvf(1,k), 1, 0.0d0, vf, 1)
        tmp = ddot(n, alpha, 1, vb(1,k), 1)
        llf = llf + log(tmp)
        nn = nn + 1.0d0
        wp(k) = 1.0d0 / tmp
        call daxpy(n, wp(k), vb(1,k), 1, eb, 1)
        call daxpy(n, wp(k), vf, 1, ey, 1)
      end if
    end do

    ! for the interval [t_m, infinity)
    if (gdatlast >= 0) then
      tmp = ddot(n, alpha, 1, barvb(1,m), 1)
      llf = llf + gdatlast * log(tmp) - logfact(gdatlast)
      nn = nn + gdatlast
      uu = uu + tmp
      wg(m+1) = gdatlast / tmp
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if

    ! compupte weights for unobserved periods
    do k = 1, m
      if (gdat(k) == -1) then
        tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
        tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)
        wg(k) = omega
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if
    end do
    if (gdatlast == -1) then
      wg(m+1) = omega
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if
    llf = llf + nn * log(omega) - omega * uu

    ! compute vectors for convolution
    vc(1:n,m) = 0.0d0
    call daxpy(n, wg(m+1)-wg(m), baralpha, 1, vc(1,m), 1)
    if (idat(m) == 1) then
      call daxpy(n, wp(m), alpha, 1, vc(1,m), 1)
    end if
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_csc_vec('T', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call daxpy(n, wg(k+1)-wg(k), baralpha, 1, vc(1,k), 1)

      if (idat(k) == 1) then
        call daxpy(n, wp(k), alpha, 1, vc(1,k), 1)
      end if
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call dcscr(n, n, wg(k+1)-wg(k), baralpha, 1, barvb(1,k), 1, &
        en, colptr, rowind, nnz)
      call mexpconv_unif_csc_vec('T', 'N', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do
    call dcscr(n, n, wg(1), baralpha, 1, barvb(1,0), 1, &
      en, colptr, rowind, nnz)

    etotal = nn + omega * (1.0d0 - uu)
    eb(1:n) = alpha(1:n) * eb(1:n)
    ey(1:n) = tau(1:n) * ey(1:n)
    do j = base, base+n-1
      do z = colptr(j), colptr(j+1)-1
        i = rowind(z)
        if (i == j) then
          ez(i) = en(z)
          exit
        end if
      end do
    end do
    en(:) = spT(:) * en(:)

    deallocate(poi)
  end subroutine gph_estep_grouppoi_unif_csc

  subroutine gph_estep_grouppoi_unif_coo(n, alpha, baralpha, tau, ev, omega, &
    spT, spP, rowind, colind, nnz, qv, poieps, &
    m, tdat, gdat, gdatlast, idat, &
    etotal, eb, ey, ez, en, llf)
    use sparse
    use spblas
    use gamma
    use poisson
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: omega
    double precision, intent(in) :: alpha(1:n), baralpha(1:n), tau(1:n), ev(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    integer, intent(in) :: rowind(base:base+nnz-1), colind(base:base+nnz-1)
    double precision, intent(in) :: qv, poieps
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, tmax
    double precision :: tmp, nn, uu
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: tildevf(1:n), tildevb(1:n)
    double precision :: barvf(1:n,0:m), barvb(1:n,0:m)
    double precision :: vf(1:n), vb(1:n,0:m)
    double precision :: vc(1:n,0:m)
    double precision :: wg(1:m+1), wp(1:m)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, baralpha, 1, barvf(1,0), 1)
    call dcopy(n, ev, 1, barvb(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_coo_vec('T', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        barvf(1,k-1), 1, barvf(1,k), 1, atol)
      call mexp_unif_coo_vec('N', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        barvb(1,k-1), 1, barvb(1,k), 1, atol)

      call dcoomv('N', n, n, -1.0d0, spT, rowind, colind, nnz, &
        barvb(1,k), 1, 0.0d0, vb(1,k), 1)

      tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
      tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)

      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmp = ddot(n, alpha, 1, tildevb, 1)
        llf = llf + gdat(k) * log(tmp) - logfact(gdat(k))
        nn = nn + gdat(k)
        uu = uu + tmp
        wg(k) = gdat(k) / tmp
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if

      if (idat(k) == 1) then
        call dcoomv('T', n, n, -1.0d0, spT, rowind, colind, nnz, &
          barvf(1,k), 1, 0.0d0, vf, 1)
        tmp = ddot(n, alpha, 1, vb(1,k), 1)
        llf = llf + log(tmp)
        nn = nn + 1.0d0
        wp(k) = 1.0d0 / tmp
        call daxpy(n, wp(k), vb(1,k), 1, eb, 1)
        call daxpy(n, wp(k), vf, 1, ey, 1)
      end if
    end do

    ! for the interval [t_m, infinity)
    if (gdatlast >= 0) then
      tmp = ddot(n, alpha, 1, barvb(1,m), 1)
      llf = llf + gdatlast * log(tmp) - logfact(gdatlast)
      nn = nn + gdatlast
      uu = uu + tmp
      wg(m+1) = gdatlast / tmp
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if

    ! compupte weights for unobserved periods
    do k = 1, m
      if (gdat(k) == -1) then
        tildevf(1:n) = barvf(1:n,k-1) - barvf(1:n,k)
        tildevb(1:n) = barvb(1:n,k-1) - barvb(1:n,k)
        wg(k) = omega
        call daxpy(n, wg(k), tildevb, 1, eb, 1)
        call daxpy(n, wg(k), tildevf, 1, ey, 1)
      end if
    end do
    if (gdatlast == -1) then
      wg(m+1) = omega
      call daxpy(n, wg(m+1), barvb(1,m), 1, eb, 1)
      call daxpy(n, wg(m+1), barvf(1,m), 1, ey, 1)
    end if
    llf = llf + nn * log(omega) - omega * uu

    ! compute vectors for convolution
    vc(1:n,m) = 0.0d0
    call daxpy(n, wg(m+1)-wg(m), baralpha, 1, vc(1,m), 1)
    if (idat(m) == 1) then
      call daxpy(n, wp(m), alpha, 1, vc(1,m), 1)
    end if
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_coo_vec('T', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call daxpy(n, wg(k+1)-wg(k), baralpha, 1, vc(1,k), 1)

      if (idat(k) == 1) then
        call daxpy(n, wp(k), alpha, 1, vc(1,k), 1)
      end if
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call dcoor(n, n, wg(k+1)-wg(k), baralpha, 1, barvb(1,k), 1, &
        en, rowind, colind, nnz)
      call mexpconv_unif_coo_vec('T', 'N', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do
    call dcoor(n, n, wg(1), baralpha, 1, barvb(1,0), 1, &
      en, rowind, colind, nnz)

    etotal = nn + omega * (1.0d0 - uu)
    eb(1:n) = alpha(1:n) * eb(1:n)
    ey(1:n) = tau(1:n) * ey(1:n)
    do z = base, base+nnz-1
      i = rowind(z)
      j = colind(z)
      if (i == j) then
        ez(i) = en(z)
      end if
    end do
    en(:) = spT(:) * en(:)

    deallocate(poi)
  end subroutine gph_estep_grouppoi_unif_coo

end module gph_estep_grouppoi

