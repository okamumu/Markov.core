!
! herlang emstep
!

module herlang_emstep
  implicit none

contains

!  Description: estep for Erlang-PH with weighted time and group/truncated data
!    alpha    (in): initial vector
!    nshape    (in): shape parameter vector
!    rate     (in): rate parameter vector
!    tdat     (in): interarrival time
!    wdat     (in): weights for interarrivals
!    gdat     (in): # of arrivals (-1 means NA)
!    gdatlast (in): # of arrivals in [lasttime, infinity] (-1 means NA)
!    idat     (in): indicator whether an arrival occurs at the last instant
!    etotal  (out): expected # of arrivals
!    eb      (out): expected # of starts
!    ew      (out): expected sojourn time?
!    return value -> llf (log-likelihood)

  subroutine herlang_estep_wtime(n, alpha, nshape, rate, &
    m, tdat, wdat, etotal, eb, ew, llf)
    use erlang_dist
    integer, intent(in) :: n, m
    double precision, intent(in) :: alpha(1:n), rate(1:n)
    integer, intent(in) :: nshape(1:n)
    double precision, intent(in) :: tdat(1:m), wdat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ew(1:n)
    double precision, intent(out) :: llf

    integer :: i, k
    double precision :: ct, scale, ddot
    double precision :: perl0(1:n,1:m), perl1(1:n,1:m)

    do k = 1, m
      ct = ct + tdat(k)
      do i = 1, n
        perl0(i,k) = pdf_erlang(nshape(i), rate(i), ct)
        perl1(i,k) = ct * perl0(i,k)
      end do
    end do

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ew = 0.0d0

    do k = 1, m
      scale = ddot(n, alpha, 1, perl0(1,k), 1)
      call daxpy(n, wdat(k)/scale, perl0(1,k), 1, eb, 1)
      call daxpy(n, wdat(k)/scale, perl1(1,k), 1, ew, 1)
      llf = llf + wdat(k) * log(scale)
      etotal = etotal + wdat(k)
    end do

    eb(1:n) = alpha(1:n) * eb(1:n)
    ew(1:n) = alpha(1:n) * ew(1:n)
  end subroutine herlang_estep_wtime

  subroutine herlang_estep_group(n, alpha, nshape, rate, &
    m, tdat, gdat, gdatlast, idat, etotal, eb, ew, llf)
    use gamma
    use erlang_dist
    integer, intent(in) :: n, m
    double precision, intent(in) :: alpha(1:n), rate(1:n)
    integer, intent(in) :: nshape(1:n)
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ew(1:n)
    double precision, intent(out) :: llf

    integer :: i, k
    double precision :: ct, scale, nn, uu, ddot
    double precision :: perl0(1:n,1:m), perl1(1:n,1:m)
    double precision :: cerl0(1:n,0:m+1), cerl1(1:n,0:m+1)
    double precision :: tmpv0(1:n), tmpv1(1:n)

    cerl0(1:n,0) = 0.0d0
    cerl1(1:n,0) = 0.0d0
    do k = 1, m
      ct = ct + tdat(k)
      do i = 1, n
        perl0(i,k) = pdf_erlang(nshape(i), rate(i), ct)
        perl1(i,k) = ct * perl0(i,k)
        cerl0(i,k) = cdf_erlang(nshape(i), rate(i), ct)
        cerl1(i,k) = (nshape(i) / rate(i)) * cdf_erlang(nshape(i)+1, rate(i), ct)
      end do
    end do
    cerl0(1:n,m+1) = 1.0d0
    cerl1(1:n,m+1) = nshape(1:n) / rate(1:n)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ew = 0.0d0
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmpv0(1:n) = cerl0(1:n,k) - cerl0(1:n,k-1)
        tmpv1(1:n) = cerl1(1:n,k) - cerl1(1:n,k-1)
        scale = ddot(n, alpha, 1, tmpv0, 1)
        nn = nn + gdat(k)
        uu = uu + scale
        call daxpy(n, gdat(k)/scale, tmpv0, 1, eb, 1)
        call daxpy(n, gdat(k)/scale, tmpv1, 1, ew, 1)
        llf = llf + gdat(k) * log(scale) - logfact(gdat(k))
      end if
      if (idat(k) == 1) then
        scale = ddot(n, alpha, 1, perl0(1,k), 1)
        nn = nn + 1.0d0
        call daxpy(n, 1.0d0/scale, perl0(1,k), 1, eb, 1)
        call daxpy(n, 1.0d0/scale, perl1(1,k), 1, ew, 1)
        llf = llf + log(scale)
      end if
    end do
    if (gdatlast >= 0) then
      tmpv0(1:n) = cerl0(1:n,m+1) - cerl0(1:n,m)
      tmpv1(1:n) = cerl1(1:n,m+1) - cerl1(1:n,m)
      scale = ddot(n, alpha, 1, tmpv0, 1)
      nn = nn + gdatlast
      uu = uu + scale
      call daxpy(n, gdatlast/scale, tmpv0, 1, eb, 1)
      call daxpy(n, gdatlast/scale, tmpv1, 1, ew, 1)
      llf = llf + gdatlast * log(scale) - logfact(gdatlast)
    end if

    do k = 1, m
      if (gdat(k) == -1) then
        tmpv0(1:n) = cerl0(1:n,k) - cerl0(1:n,k-1)
        tmpv1(1:n) = cerl1(1:n,k) - cerl1(1:n,k-1)
        call daxpy(n, nn/uu, tmpv0, 1, eb, 1)
        call daxpy(n, nn/uu, tmpv1, 1, ew, 1)
      end if
    end do
    if (gdatlast == -1) then
      tmpv0(1:n) = cerl0(1:n,m+1) - cerl0(1:n,m)
      tmpv1(1:n) = cerl1(1:n,m+1) - cerl1(1:n,m)
      call daxpy(n, nn/uu, tmpv0, 1, eb, 1)
      call daxpy(n, nn/uu, tmpv1, 1, ew, 1)
    end if
    llf = llf + loggamma(nn + 1.0d0) - nn * log(uu)

    etotal = nn/uu
    eb(1:n) = alpha(1:n) * eb(1:n)
    ew(1:n) = alpha(1:n) * ew(1:n)
  end subroutine herlang_estep_group

  subroutine herlang_estep_grouppoi(n, alpha, nshape, rate, omega, &
    m, tdat, gdat, gdatlast, idat, etotal, eb, ew, llf)
    use gamma
    use erlang_dist
    integer, intent(in) :: n, m
    double precision, intent(in) :: alpha(1:n), rate(1:n), omega
    integer, intent(in) :: nshape(1:n)
    double precision, intent(in) :: tdat(1:m)
    integer, intent(in) :: gdat(1:m), gdatlast, idat(1:m)
    double precision, intent(out) :: etotal, eb(1:n), ew(1:n)
    double precision, intent(out) :: llf

    integer :: i, k
    double precision :: ct, scale, nn, uu, ddot
    double precision :: perl0(1:n,1:m), perl1(1:n,1:m)
    double precision :: cerl0(1:n,0:m+1), cerl1(1:n,0:m+1)
    double precision :: tmpv0(1:n), tmpv1(1:n)

    cerl0(1:n,0) = 0.0d0
    cerl1(1:n,0) = 0.0d0
    do k = 1, m
      ct = ct + tdat(k)
      do i = 1, n
        perl0(i,k) = pdf_erlang(nshape(i), rate(i), ct)
        perl1(i,k) = ct * perl0(i,k)
        cerl0(i,k) = cdf_erlang(nshape(i), rate(i), ct)
        cerl1(i,k) = (nshape(i) / rate(i)) * cdf_erlang(nshape(i)+1, rate(i), ct)
      end do
    end do
    cerl0(1:n,m+1) = 1.0d0
    cerl1(1:n,m+1) = nshape(1:n) / rate(1:n)

    llf = 0.0d0
    etotal = 0.0d0
    eb = 0.0d0
    ew = 0.0d0
    nn = 0.0d0
    uu = 0.0d0

    do k = 1, m
      if (gdat(k) >= 0 .and. tdat(k) /= 0.0d0) then
        tmpv0(1:n) = cerl0(1:n,k) - cerl0(1:n,k-1)
        tmpv1(1:n) = cerl1(1:n,k) - cerl1(1:n,k-1)
        scale = ddot(n, alpha, 1, tmpv0, 1)
        nn = nn + gdat(k)
        uu = uu + scale
        call daxpy(n, gdat(k)/scale, tmpv0, 1, eb, 1)
        call daxpy(n, gdat(k)/scale, tmpv1, 1, ew, 1)
        llf = llf + gdat(k) * log(scale) - logfact(gdat(k))
      end if
      if (idat(k) == 1) then
        scale = ddot(n, alpha, 1, perl0(1,k), 1)
        nn = nn + 1.0d0
        call daxpy(n, 1.0d0/scale, perl0(1,k), 1, eb, 1)
        call daxpy(n, 1.0d0/scale, perl1(1,k), 1, ew, 1)
        llf = llf + log(scale)
      end if
    end do
    if (gdatlast >= 0) then
      tmpv0(1:n) = cerl0(1:n,m+1) - cerl0(1:n,m)
      tmpv1(1:n) = cerl1(1:n,m+1) - cerl1(1:n,m)
      scale = ddot(n, alpha, 1, tmpv0, 1)
      nn = nn + gdatlast
      uu = uu + scale
      call daxpy(n, gdatlast/scale, tmpv0, 1, eb, 1)
      call daxpy(n, gdatlast/scale, tmpv1, 1, ew, 1)
      llf = llf + gdatlast * log(scale) - logfact(gdatlast)
    end if

    do k = 1, m
      if (gdat(k) == -1) then
        tmpv0(1:n) = cerl0(1:n,k) - cerl0(1:n,k-1)
        tmpv1(1:n) = cerl1(1:n,k) - cerl1(1:n,k-1)
        call daxpy(n, omega, tmpv0, 1, eb, 1)
        call daxpy(n, omega, tmpv1, 1, ew, 1)
      end if
    end do
    if (gdatlast == -1) then
      tmpv0(1:n) = cerl0(1:n,m+1) - cerl0(1:n,m)
      tmpv1(1:n) = cerl1(1:n,m+1) - cerl1(1:n,m)
      call daxpy(n, omega, tmpv0, 1, eb, 1)
      call daxpy(n, omega, tmpv1, 1, ew, 1)
    end if
    llf = llf + nn * log(omega) - omega * uu

    etotal = nn + omega * (1.0d0 - uu)
    eb(1:n) = alpha(1:n) * eb(1:n)
    ew(1:n) = alpha(1:n) * ew(1:n)
  end subroutine herlang_estep_grouppoi

  subroutine herlang_mstep(n, etotal, eb, ew, alpha, nshape, rate)
    integer, intent(in) :: n
    integer, intent(in) :: nshape(1:n)
    double precision, intent(in) :: etotal, eb(1:n), ew(1:n)
    double precision, intent(out) :: alpha(1:n), rate(1:n)

    alpha(1:n) = eb(1:n) / etotal
    rate(1:n) = nshape(1:n) * eb(1:n) / ew(1:n)

  end subroutine herlang_mstep

end module herlang_emstep

