!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ju@Sussex
! FORTRAN routines to compute low level products while generating
! the scanning strategy.
! Main purposes is interfacing with python (using f2py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module scanning_strategy_f

contains

    subroutine mapmaking(pix_global, nhit_loc, npix, num_pts)
        implicit none
        ! Simple map-making: project time ordered data into sky maps for
        ! visualisation.

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: num_pts, npix
        integer(I4B), intent(in) :: pix_global(0: num_pts - 1)
        integer(I4B), intent(inout) :: nhit_loc(0: npix - 1)

        ! LOCAL
        integer(I4B)             :: i, pix

        do i=0, num_pts - 1
            pix = pix_global(i)

            ! Number of hits per pixel
            nhit_loc(pix) = nhit_loc(pix) + 1
        enddo

    end subroutine

    subroutine convolve_focalplane_f(bore_nhits, focalplane_nhits,&
    pixels, bolo_per_pix, npix_loc, boost)
        implicit none
        ! Given a number of hits map perform the focal plane convolution.
        ! Original author (python routine): Neil Goeckner-Wald.
        ! Modifications by Julien Peloton (including fortran interface).

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: npix_loc
        integer(I4B), intent(in) :: pixels(0: npix_loc - 1)
        real(DP), intent(in)     :: bore_nhits
        real(DP), intent(in)     :: bolo_per_pix, boost
        real(DP), intent(inout)  :: focalplane_nhits(0: npix_loc - 1)

        ! LOCAL
        integer(I4B)             :: i, pix

        ! Loop
        do i=0, npix_loc - 1
            pix = pixels(i)
            focalplane_nhits(pix) = focalplane_nhits(pix) + &
            bore_nhits * bolo_per_pix * boost
        enddo

    end subroutine

    subroutine run_one_scan_f(pb_az_array, pb_mjd_array, running_az, &
        upper_az, lower_az, az_speed, pb_az_dir, second, sampling_freq, num_pts)
        implicit none
        ! Generate one observation (i.e. one CES) of the telescope.

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: num_pts
        real(DP), intent(in)     :: upper_az, lower_az, az_speed
        real(DP), intent(in)     :: second, sampling_freq
        real(DP), intent(inout)  :: running_az, pb_az_dir
        real(DP), intent(inout)  :: pb_az_array(0: num_pts - 1)
        real(DP), intent(inout)  :: pb_mjd_array(0: num_pts - 1)

        ! LOCAL
        integer(I4B)             :: t

        ! Loop. Start at one since the initialisation
        ! of the scan is done outside.
        do t=1, num_pts - 1
            ! Set the Azimuth and time
            pb_az_array(t) = running_az

            ! Case to change the direction of the scan
            if (running_az .gt. upper_az) then
                pb_az_dir = -1.0
            elseif (running_az .lt. lower_az) then
                pb_az_dir = 1.0
            endif

            running_az = running_az + az_speed * pb_az_dir / sampling_freq

            ! Increment the time by one second / sampling rate
            pb_mjd_array(t) = pb_mjd_array(t - 1) + second / sampling_freq
        end do

    end subroutine
end module
