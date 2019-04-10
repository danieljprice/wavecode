!
!--routines to return eta/zeta for a spinning black hole
!
module blackhole
 implicit none
 real    :: a_spin,rs
 logical :: bhset = .false.
 integer :: ifreq

 integer, parameter :: ifreq_gr       = 1
 integer, parameter :: ifreq_grapprox = 2
 integer, parameter :: ifreq_pn       = 3

contains

subroutine set_bh(spin,rsch)
 real, intent(in) :: spin,rsch

 a_spin = spin
 rs     = rsch
 bhset  = .true.

end subroutine set_bh

subroutine get_bh(r,eta,zeta)
 real,    intent(in)  :: r
 real,    intent(out) :: eta,zeta
 real :: term

 if (.not. bhset) stop 'error: black hole not set'
 term = 0.

 select case(ifreq)

 !-- GR frequencies
 case(ifreq_gr)
    eta  = -1.5*rs/r + sqrt(2.)*a_spin*(rs/r)**1.5 - 3./8.*(a_spin*rs/r)**2
    zeta = -a_spin/sqrt(2.)*sqrt((rs/r)**3) + 3./8.*(a_spin*rs/r)**2

 !-- Approximate GR frequencies (as in LOP02)
 case(ifreq_grapprox)
    eta = -1.5*rs/r
    zeta = -a_spin/sqrt(2.)*sqrt((rs/r)**3)

 !-- Post-Newtonian frequencies
 case(ifreq_pn)
    ! In the PN case, we assume G=M=c=1. See frequencies in NPN15.
    print*,'Using the PN frequencies!'
    term = 2.*r**1.5 - 4.*a_spin
    eta  = 3.*a_spin/term
    zeta = -4.*a_spin/term
 case default
    STOP 'Bad choice of ifreq'
 end select

end subroutine get_bh

end module blackhole
