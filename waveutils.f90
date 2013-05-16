module waveutils
 implicit none
 integer, parameter :: maxgrid = 4004
 complex*16, dimension(maxgrid) :: za1,za2,zd1,zd2  ! zvbles
 complex*16                     :: zi               ! zvbles
 real*8, dimension(2*maxgrid) :: r,dr,rsq,r12,r32 ! grid
 real*8, dimension(2*maxgrid) :: omega,rho,csq    ! disc
 real*8, dimension(2*maxgrid) :: eta, zeta        ! precess
 real*8  :: rin,rout           ! grid
 real*8  :: alpha           ! diss
 real*8  :: rstep,wstep     ! ics
 real*8  :: time,dt,ctime   ! tempus
 integer :: nstep,nfile     ! tempus
 integer :: n
 real*8, parameter :: pi = 4.*atan(1.)
 real*8  :: etazero,zetazero,honr         ! consts
 character(len=20) :: mode

end module waveutils
