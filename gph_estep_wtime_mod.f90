!
! gph estep
!

module gph_estep_wtime
  implicit none

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

  subroutine gph_estep_wtime_unif_dense(n, alpha, tau, &
    T, P, qv, poieps, m, tdat, wdat, etotal, eb, ey, ez, en, llf)
    use poisson
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, m
    double precision, intent(in) :: alpha(1:n), tau(1:n)
    double precision, intent(in) :: T(1:n,1:n), P(1:n,1:n), qv, poieps
    double precision, intent(in) :: tdat(1:m), wdat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(1:n)
    double precision, intent(out) :: en(1:n,1:n)
    double precision, intent(out) :: llf

    integer :: k, right
    double precision :: weight, scale, tllf, tmax
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: vf(1:n,0:m), vb(1:n,0:m), vc(1:n,0:m)
    double precision :: blf(1:m)

    llf = 0.0d0
    tllf = 0.0d0

    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, alpha, 1, vf(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_dense_vec('T', n, P, n, qv, left, right, poi, weight, &
        vf(1,k-1), 1, vf(1,k), 1, atol)
      scale = ddot(n, vf(1,k), 1, tau, 1)
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
      call daxpy(n, wdat(k), vf(1,k), 1, ey, 1)
      blf(k) = scale

      call mexp_unif_dense_vec('N', n, P, n, qv, left, right, poi, weight, &
        vb(1,k-1), 1, vb(1,k), 1, atol)
      scale = ddot(n, alpha, 1, vb(1,k), 1)
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      call daxpy(n, wdat(k), vb(1,k), 1, eb, 1)

      etotal = etotal + wdat(k)
      tllf = tllf + log(blf(k))
      llf = llf + wdat(k) * tllf
    end do

    vc(1:n,m) = 0.0d0
    call daxpy(n, wdat(m)/blf(m), alpha, 1, vc(1,m), 1)
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_dense_vec('T', n, P, n, qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call dscal(n, 1.0/blf(k), vc(1,k), 1)
      call daxpy(n, wdat(k)/blf(k), alpha, 1, vc(1,k), 1)
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexpconv_unif_dense_vec('T', 'N', n, P, n, qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, n, atol)
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    ey(1:n) = tau(1:n) * ey(1:n)
    call dcopy(n, en(1,1), n+1, ez, 1)
    en(1:n,1:n) = T(1:n,1:n) * en(1:n,1:n)

    deallocate(poi)
  end subroutine gph_estep_wtime_unif_dense

  subroutine gph_estep_wtime_unif_csr(n, alpha, tau, &
    spT, spP, rowptr, colind, nnz, qv, poieps, &
    m, tdat, wdat, etotal, eb, ey, ez, en, llf)
    use poisson
    use sparse
    use spblas
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: alpha(1:n), tau(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    double precision, intent(in) :: qv, poieps
    integer, intent(in) :: rowptr(base:base+n), colind(base:base+nnz-1)
    double precision, intent(in) :: tdat(1:m), wdat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, scale, tllf, tmax
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: vf(1:n,0:m), vb(1:n,0:m), vc(1:n,0:m)
    double precision :: blf(1:m)

    llf = 0.0d0
    tllf = 0.0d0

    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, alpha, 1, vf(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_csr_vec('T', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vf(1,k-1), 1, vf(1,k), 1, atol)
      scale = ddot(n, vf(1,k), 1, tau, 1)
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
      call daxpy(n, wdat(k), vf(1,k), 1, ey, 1)
      blf(k) = scale

      call mexp_unif_csr_vec('N', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vb(1,k-1), 1, vb(1,k), 1, atol)
      scale = ddot(n, alpha, 1, vb(1,k), 1)
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      call daxpy(n, wdat(k), vb(1,k), 1, eb, 1)

      etotal = etotal + wdat(k)
      tllf = tllf + log(blf(k))
      llf = llf + wdat(k) * tllf
    end do

    vc(1:n,m) = 0.0d0
    call daxpy(n, wdat(m)/blf(m), alpha, 1, vc(1,m), 1)
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_csr_vec('T', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call dscal(n, 1.0/blf(k), vc(1,k), 1)
      call daxpy(n, wdat(k)/blf(k), alpha, 1, vc(1,k), 1)
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexpconv_unif_csr_vec('T', 'N', n, spP, rowptr, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do

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
  end subroutine gph_estep_wtime_unif_csr

  subroutine gph_estep_wtime_unif_csc(n, alpha, tau, &
    spT, spP, colptr, rowind, nnz, qv, poieps, &
    m, tdat, wdat, etotal, eb, ey, ez, en, llf)
    use poisson
    use sparse
    use spblas
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: alpha(1:n), tau(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    double precision, intent(in) :: qv, poieps
    integer, intent(in) :: colptr(base:base+n), rowind(base:base+nnz-1)
    double precision, intent(in) :: tdat(1:m), wdat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, scale, tllf, tmax
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: vf(1:n,0:m), vb(1:n,0:m), vc(1:n,0:m)
    double precision :: blf(1:m)

    llf = 0.0d0
    tllf = 0.0d0

    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, alpha, 1, vf(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_csc_vec('T', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vf(1,k-1), 1, vf(1,k), 1, atol)
      scale = ddot(n, vf(1,k), 1, tau, 1)
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
      call daxpy(n, wdat(k), vf(1,k), 1, ey, 1)
      blf(k) = scale

      call mexp_unif_csc_vec('N', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vb(1,k-1), 1, vb(1,k), 1, atol)
      scale = ddot(n, alpha, 1, vb(1,k), 1)
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      call daxpy(n, wdat(k), vb(1,k), 1, eb, 1)

      etotal = etotal + wdat(k)
      tllf = tllf + log(blf(k))
      llf = llf + wdat(k) * tllf
    end do

    vc(1:n,m) = 0.0d0
    call daxpy(n, wdat(m)/blf(m), alpha, 1, vc(1,m), 1)
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_csc_vec('T', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call dscal(n, 1.0/blf(k), vc(1,k), 1)
      call daxpy(n, wdat(k)/blf(k), alpha, 1, vc(1,k), 1)
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexpconv_unif_csc_vec('T', 'N', n, spP, colptr, rowind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do

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
  end subroutine gph_estep_wtime_unif_csc

  subroutine gph_estep_wtime_unif_coo(n, alpha, tau, &
    spT, spP, rowind, colind, nnz, qv, poieps, &
    m, tdat, wdat, etotal, eb, ey, ez, en, llf)
    use poisson
    use sparse
    use spblas
    use mexp_unif
    use mexpconv_unif
    integer, parameter :: base = sparse_base_index
    integer, parameter :: left = 0
    double precision, parameter :: atol = 0.0d0

    integer, intent(in) :: n, nnz, m
    double precision, intent(in) :: alpha(1:n), tau(1:n)
    double precision, intent(in) :: spT(1:nnz), spP(1:nnz)
    double precision, intent(in) :: qv, poieps
    integer, intent(in) :: rowind(base:base+nnz-1), colind(base:base+nnz-1)
    double precision, intent(in) :: tdat(1:m), wdat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ey(1:n), ez(base:base+n-1)
    double precision, intent(out) :: en(base:base+nnz-1)
    double precision, intent(out) :: llf

    integer :: i, j, k, z, right
    double precision :: weight, scale, tllf, tmax
    double precision :: ddot
    double precision, allocatable :: poi(:)
    double precision :: vf(1:n,0:m), vb(1:n,0:m), vc(1:n,0:m)
    double precision :: blf(1:m)

    llf = 0.0d0
    tllf = 0.0d0

    etotal = 0.0d0
    eb = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    en = 0.0d0

    tmax = maxval(tdat(1:m))
    right = poisson_rightbound(qv*tmax, poieps) + 1
    allocate(poi(0:right))

    call dcopy(n, alpha, 1, vf(1,0), 1)
    call dcopy(n, tau, 1, vb(1,0), 1)
    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexp_unif_coo_vec('T', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vf(1,k-1), 1, vf(1,k), 1, atol)
      scale = ddot(n, vf(1,k), 1, tau, 1)
      call dscal(n, 1.0d0/scale, vf(1,k), 1)
      call daxpy(n, wdat(k), vf(1,k), 1, ey, 1)
      blf(k) = scale

      call mexp_unif_coo_vec('N', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vb(1,k-1), 1, vb(1,k), 1, atol)
      scale = ddot(n, alpha, 1, vb(1,k), 1)
      call dscal(n, 1.0d0/scale, vb(1,k), 1)
      call daxpy(n, wdat(k), vb(1,k), 1, eb, 1)

      etotal = etotal + wdat(k)
      tllf = tllf + log(blf(k))
      llf = llf + wdat(k) * tllf
    end do

    vc(1:n,m) = 0.0d0
    call daxpy(n, wdat(m)/blf(m), alpha, 1, vc(1,m), 1)
    do k = m-1, 1, -1
      right = poisson_rightbound(qv*tdat(k+1), poieps) + 1
      call poisson_prob(qv*tdat(k+1), left, right, poi, weight)

      call mexp_unif_coo_vec('T', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k+1), 1, vc(1,k), 1, atol)
      call dscal(n, 1.0/blf(k), vc(1,k), 1)
      call daxpy(n, wdat(k)/blf(k), alpha, 1, vc(1,k), 1)
    end do

    do k = 1, m
      right = poisson_rightbound(qv*tdat(k), poieps) + 1
      call poisson_prob(qv*tdat(k), left, right, poi, weight)

      call mexpconv_unif_coo_vec('T', 'N', n, spP, rowind, colind, nnz, &
        qv, left, right, poi, weight, &
        vc(1,k), 1, vb(1,k-1), 1, vb(1,k-1), 1, en, atol)
    end do

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
  end subroutine gph_estep_wtime_unif_coo

end module gph_estep_wtime

