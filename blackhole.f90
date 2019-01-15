!
!--routines to return eta/zeta for a spinning black hole
!
module blackhole
 implicit none
 real :: a_spin,rs
 logical :: bhset = .false.

contains

 subroutine set_bh(spin,rsch)
  real, intent(in) :: spin,rsch

  a_spin = spin
  rs = rsch
  bhset = .true.

 end subroutine set_bh

 subroutine get_bh(r,eta,zeta)
  real, intent(in)  :: r
  real, intent(out) :: eta,zeta

  if (.not. bhset) stop 'error: black hole not set'

  eta = -1.5*rs/r
  zeta = -a_spin/sqrt(2.)*sqrt((rs/r)**3)

 end subroutine get_bh

end module blackhole
