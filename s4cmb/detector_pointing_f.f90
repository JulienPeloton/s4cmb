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
        integer(I4B)             :: angle

        do angle=0, 4 * n - 1, 4

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

    subroutine quat_to_radecpa_fortran_f(q0, q1, q2, q3, phi, theta, psi, n)
        implicit none
        ! Routine to compute phi/theta/psi from a sequence
        ! of quaternions. Computation is done in fortran.
        !
        ! WARNING: you still need to convert phi/theta/psi to get to RA/Dec/PA.
        !
        ! Parameters
        ! ----------
        ! seq : array of arrays
        !     Array of quaternions.
        !
        ! Returns
        ! ----------
        ! phi : 1d array
        ! theta : 1d array
        ! psi : 1d array

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        ! F2PY params
        integer(I4B), intent(in) :: n
        real(DP), intent(in)     :: q0(0 : n - 1), q1(0 : n - 1)
        real(DP), intent(in)     :: q2(0 : n - 1), q3(0 : n - 1)
        real(DP), intent(inout)  :: phi(0 : n - 1), theta(0 : n - 1), psi(0 : n - 1)

        ! LOCAL
        integer(I4B)             :: i

        do i=0, n - 1
            phi(i) = atan2(2.0 * (q0(i) * q1(i) + q2(i) * q3(i)), &
                1.0 - 2.0 * (q1(i) * q1(i) + q2(i) * q2(i)))
            theta(i) = asin(2.0 * (q0(i) * q2(i) - q3(i) * q1(i)))
            psi(i) = atan2(2.0 * (q0(i) * q3(i) + q1(i) * q2(i)), &
                1.0 - 2.0 * (q2(i) * q2(i) + q3(i) * q3(i)))
        enddo

    end subroutine

end module
