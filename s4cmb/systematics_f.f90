!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ju@Sussex
! FORTRAN routines to compute low level products while injecting systematics.
! Main purposes is interfacing with python (using f2py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module systematics_f

contains

    subroutine inject_crosstalk_inside_squid_f(tsout, local_indices, global_indices, &
        radius, cross_amp, beta, num_chans_local, num_chans_global, nts)
        implicit none
        ! Introduce leakage between neighboring bolometers within a SQUID.
        ! local_indices is an array with bolometers indices within their SQUID.
        ! global_indices is an array with bolometers indices within the focal plane.
        ! For other parameters and docs, see the python routine (s4cmb)
        ! systematics.inject_crosstalk_inside_squid

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: num_chans_local, num_chans_global
        integer(I4B), intent(in) :: beta, radius, nts
        integer(I4B), intent(in) :: local_indices(0: num_chans_local - 1)
        integer(I4B), intent(in) :: global_indices(0: num_chans_local - 1)
        real(DP),     intent(in) :: cross_amp(0: num_chans_global - 1)
        real(DP),  intent(inout) :: tsout(0 : num_chans_global - 1, 0 : nts - 1)

        ! LOCAL
        integer(I4B)             :: local_index, local_index2
        integer(I4B)             :: global_index, global_index2
        integer(I4B)             :: i, i2, separation_length, n

        ! Pick up a channel...
        do i=0, num_chans_local - 1
            local_index = local_indices(i)
            global_index = global_indices(i)

            ! ... And oop over all the others.
            do i2=0, num_chans_local - 1
                local_index2 = local_indices(i2)
                global_index2 = global_indices(i2)
                separation_length = ABS(local_index - local_index2)

                ! Look whether channels can cross talk.
                if (separation_length .gt. 0 .and. separation_length .lt. radius + 1) then

                    ! Here is your model of crosstalk.
                    tsout(global_index, :) = tsout(global_index, :) + &
                        cross_amp(global_index2) / &
                        separation_length**beta * tsout(global_index2, :)
                endif
            enddo
        enddo

    end subroutine

end module
