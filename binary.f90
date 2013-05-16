!
!--routines to return eta/zeta for a binary potential
!
module binary
 implicit none
 real*8 :: m1,m2,rs1,rs2
 logical :: binaryset = .false.
 
contains

 subroutine set_binary(mass1,mass2,r1,r2)
  real*8, intent(in) :: mass1,mass2,r1,r2
  
  m1 = mass1
  m2 = mass2
  rs1 = r1
  rs2 = r2
  binaryset = .true.

 end subroutine set_binary

 subroutine get_binary(r,eta,zeta,omega)
  real*8, intent(in)  :: r
  real*8, intent(out) :: eta,zeta,omega
  real*8 :: term1,term2,omega2,omegaz2,kappa2
  
  if (.not. binaryset) stop 'error: binary not set'

  term1 = (m1 + m2)/r**3
  term2 = (m1*rs1**2 + m2*rs2**2)/r**5
  omega2  = term1 +    0.75*term2
  omegaz2 = term1 + 3.*0.75*term2
  kappa2  = term1 -    0.75*term2

  eta  = (kappa2 - omega2)/(2.*omega2)
  zeta = (omegaz2 - omega2)/(2.*omega2)
  omega = sqrt(omega2)

 end subroutine get_binary
 
end module binary
