!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ju@Sussex
! FORTRAN routines to compute low level products
! Main purposes is interfacing with python (using f2py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module tod_f

contains

    subroutine tod2map_alldet_f(d, w, dc, ds, cc, cs, ss, nhit, waferi1d, &
    waferpa, waferts, diff_weight, sum_weight, npix, nt, &
    wafermask_pixel, nskypix)
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: npix, nt, nskypix
        integer(I4B), intent(in) :: waferi1d(0:npix*nt - 1)
        integer(I4B), intent(in) :: wafermask_pixel(0:npix*nt - 1)
        real(DP), intent(in)     :: waferpa(0:npix*nt - 1), waferts(0:npix*nt*2 - 1)
        real(DP), intent(in)     :: diff_weight(0:npix - 1), sum_weight(0:npix - 1)

        real(DP), intent(inout)  :: d(0:nskypix - 1), w(0:nskypix - 1), dc(0:nskypix - 1)
        real(DP), intent(inout)  :: ds(0:nskypix - 1), cc(0:nskypix - 1)
        real(DP), intent(inout)  :: cs(0:nskypix - 1), ss(0:nskypix - 1)
        integer(I4B), intent(inout) :: nhit(0:nskypix - 1)

        integer(I4B)             :: i, j, ipix, pixel
        integer(I4B)             :: ict, icb
        real(DP)                 :: sum, diff, c, s

        do j=0, npix - 1
            do i=0, nt - 1
                ipix = i + j * nt
                if (wafermask_pixel(ipix) .gt. 0 .and. waferi1d(ipix) .gt. 0) then
                    ict = i + 2*j*nt
                    icb = i + (2*j + 1)*nt

                    pixel = waferi1d(ipix)

                    sum = 0.5*(waferts(ict) + waferts(icb))
                    diff = 0.5*(waferts(ict) - waferts(icb))
                    c = cos(2.0*waferpa(ipix))
                    s = sin(2.0*waferpa(ipix))

                    nhit(pixel) = nhit(pixel) + 1
                    w(pixel) = w(pixel) + sum_weight(j)
                    d(pixel) = d(pixel) + sum * sum_weight(j)

                    dc(pixel) = dc(pixel) + c * diff * diff_weight(j)
                    ds(pixel) = ds(pixel) + s * diff * diff_weight(j)
                    cc(pixel) = cc(pixel) + c * c * diff_weight(j)
                    cs(pixel) = cs(pixel) + c * s * diff_weight(j)
                    ss(pixel) = ss(pixel) + s * s * diff_weight(j)
                endif
            enddo
        enddo

    end subroutine

    subroutine polarized_coadd_hwp_f(d0, d4r, d4i, w0, w4, nhit, waferi1d, &
    waferpa, waferts, weight4, weight0, nch, nt, &
    wafermask_pixel, nts, nces, nskypix)
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: nch, nt, nces, nskypix
        integer(I4B), intent(in) :: waferi1d(0:nch*nt - 1), wafermask_pixel(0:nch*nt - 1)
        integer(I4B), intent(in) :: nts(0:nces - 1)
        real(DP), intent(in)     :: waferpa(0:nch*nt - 1), waferts(0:nch*nt*3 - 1)
        real(DP), intent(in)     :: weight0(0:nch - 1), weight4(0:nch - 1)

        real(DP), intent(inout)  :: d0(0:nskypix - 1), d4r(0:nskypix - 1), d4i(0:nskypix - 1)
        real(DP), intent(inout)  :: w0(0:nskypix - 1), w4(0:nskypix - 1)
        integer(I4B), intent(inout) :: nhit(0:nskypix - 1)

        integer(I4B)             :: i, j, ipix, ic, iw
        integer(I4B)             :: pixel, istop
        integer(I4B)             :: if0, i4r, i4i
        real(DP)                 :: c, s

        do j=0, nch - 1
            istop=0
            do ic=0, nces - 1
                iw = ic + j*nces
                do i=istop, nts(ic)+istop - 1
                    ipix = i + j*nt
                    if (wafermask_pixel(ipix) .gt. 0 .and. waferi1d(ipix) .gt. 0) then
                        if0 = i + j*3*nt
                        i4r = i + nt + j*3*nt
                        i4i = i + nt*2 + j*3*nt

                        pixel = waferi1d(ipix)

                        c = cos(2.0*waferpa(ipix))
                        s = sin(2.0*waferpa(ipix))

                        nhit(pixel) = nhit(pixel) + 1

                        w0(pixel) = w0(pixel) + weight0(iw)
                        w4(pixel) = w4(pixel) + weight4(iw)
                        d0(pixel) = d0(pixel)+ waferts(if0) * weight0(iw)
                        d4r(pixel) = d4r(pixel) + (c*waferts(i4r)+s*waferts(i4i)) * weight4(iw)
                        d4i(pixel) = d4i(pixel) + (s*waferts(i4r)-c*waferts(i4i)) * weight4(iw)
                    endif
                enddo
                istop = istop + nts(ic)
            enddo
        enddo
    end subroutine

end module
