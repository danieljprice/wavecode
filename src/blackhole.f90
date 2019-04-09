!
!--routines to return eta/zeta for a spinning black hole
!
module blackhole
 implicit none
 real    :: a_spin,rs
 logical :: bhset = .false.

contains

subroutine set_bh(spin,rsch)
 real, intent(in) :: spin,rsch

 a_spin = spin
 rs     = rsch
 bhset  = .true.

end subroutine set_bh

subroutine get_bh(r,eta,zeta,use_pn)
 real,    intent(in)  :: r
 real,    intent(out) :: eta,zeta
 logical, intent(in)  :: use_pn ! Set to true to use the PN frequencies rather than the GR
 real :: term

 if (.not. bhset) stop 'error: black hole not set'
 term = 0.

 if (.not.use_pn) then
    eta  = -1.5*rs/r + sqrt(2.)*a_spin*(rs/r)**1.5 - 3./8.*(a_spin*rs/r)**2
    zeta = -a_spin/sqrt(2.)*sqrt((rs/r)**3) + 3./8.*(a_spin*rs/r)**2
 else
    ! In the PN case, we assume G=M=c=1. See frequencies in NPN15.
    print*,'Using the PN frequencies!'
    term = 2.*r**1.5 - 4.*a_spin
    eta  = 3.*a_spin/term
    zeta = -4.*a_spin/term
 endif

end subroutine get_bh

end module blackhole
