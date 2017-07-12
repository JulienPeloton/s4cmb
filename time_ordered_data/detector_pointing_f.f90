!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ju@Sussex
! FORTRAN routines to compute low level products while computing detector
! pointing.
! Main purposes is interfacing with python (using f2py)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module detector_pointing_f

contains

    subroutine mult_fortran_f(p, q, pq, n)
        implicit none
        ! Multiply arrays of quaternions, when p is an array of quaternions
        ! and q is a single quaternion.

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: n
        real(DP), intent(in)     :: p(0 : 4 * n - 1)
        real(DP), intent(in)     :: q(0 : 3)
        real(DP), intent(inout)  :: pq(0 : 4 * n - 1)

        ! LOCAL
        integer(I4B)             :: pos, angle

        do angle=0, 4 * n - 1, 4
            ! angle = 4 * pos

            pq(angle + 3) = p(angle + 3) * q(3)
            pq(angle + 3) = pq(angle + 3) - (p(angle) * q(0) + &
                p(angle + 1) * q(1) + p(angle + 2) * q(2))

            pq(angle) = p(angle + 3) * q(0) + p(angle) * q(3) + &
                p(angle + 1) * q(2) - p(angle + 2) * q(1)

            pq(angle + 1) = p(angle + 3) * q(1) + p(angle + 1) * q(3)  +  &
                p(angle + 2) * q(0) - p(angle + 0) * q(2)

            pq(angle + 2) = p(angle + 3) * q(2) + p(angle + 2) * q(3)  +  &
                p(angle + 0) * q(1) - p(angle + 1) * q(0)

        enddo

    end subroutine

end module
