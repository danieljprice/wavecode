module waveutils
 implicit none
 integer, parameter :: maxgrid = 4004
 complex, dimension(maxgrid) :: za1,za2,zd1,zd2       ! zvbles
 complex                     :: zi                    ! zvbles
 real, dimension(2*maxgrid) :: r,dr,rsq,r12,r32       ! grid
 real, dimension(2*maxgrid) :: sigma,scale_height     ! more grid
 real, dimension(2*maxgrid) :: omega,rho,csq,alpha    ! disc
 real, dimension(2*maxgrid) :: eta, zeta              ! precess

 real  :: rstep,wstep       ! ics
 real  :: time,dt,ctime     ! tempus
 integer :: nstep,nfile     ! tempus

 real, parameter :: pi = 4.*atan(1.)
 real  :: etazero,zetazero         ! consts
 logical :: use_ext_sigma_profile

 !--Runtime parameters
 integer :: n
 real    :: rin,rout,alphaSS  ! grid
 real    :: honr  ! const
 character(len=20) :: mode

end module waveutils
