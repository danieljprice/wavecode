!
!  this is prog08.f
!    .... with the precession put in by leapfrog 
!
!     this solves for the propagation of a warp wave through a Keplerian
!     disc. It uses the warp equations as given in Lubow, Ogilvie &
!     Pringle (2002) MNRAS 337, 706, Section 4, using the variables A and
!     D. For a Keplerian disc we set 'zeta' and 'eta' to be zero.
!
!     for eta and zeta zero, the precession via leapfrog is not
!     required. this means that much of the complexity can be removed.
!     we no longer need both za1 and za1, nor both zd1 and zd2.
!
!  we can also add dissipation in the form of a decay of 'A' the radial
!  velocity component, using 'alpha'
!
!  Original code by Jim Pringle
!  Modernised, converted to F90 etc. by Daniel Price, April 2013
!
module waveutils
 implicit none
 integer, parameter :: maxgrid = 2002
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

end module waveutils

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

 subroutine get_binary(r,eta,zeta)
  real*8, intent(in)  :: r
  real*8, intent(out) :: eta,zeta
  real*8 :: term1,term2,omega2,omegaz2,kappa2
  
  if (.not. binaryset) stop 'error: binary not set'

  term1 = (m1 + m2)/r**3
  term2 = (m1*rs1**2 + m2*rs2**2)/r**5
  omega2  = term1 +    0.75*term2
  omegaz2 = term1 + 3.*0.75*term2
  kappa2  = term1 -    0.75*term2

  eta  = (kappa2 - omega2)/(2.*omega2)
  zeta = (omegaz2 - omega2)/(2.*omega2)

 end subroutine get_binary

end module binary

!
!--routines to return eta/zeta for a spinning black hole
!
module blackhole
 implicit none
 real*8 :: a_spin,rs
 logical :: bhset = .false.
 
contains

 subroutine set_bh(spin,rsch)
  real*8, intent(in) :: spin,rsch
  
  a_spin = spin
  rs = rsch
  bhset = .true.

 end subroutine set_bh

 subroutine get_bh(r,eta,zeta)
  real*8, intent(in)  :: r
  real*8, intent(out) :: eta,zeta
  
  if (.not. bhset) stop 'error: black hole not set'

  eta = -1.5*rs/r
  zeta = -a_spin/sqrt(2.)*sqrt((rs/r)**3)

 end subroutine get_bh

end module blackhole

!----------------
!  MAIN CODE
!----------------
program wave
 use waveutils
 use binary,    only:set_binary,get_binary
 use blackhole, only:set_bh,get_bh
 implicit none
 character(len=20) :: mode
 real*8     :: tcheck,tprint,tstop
 real*8     :: toutfile,tcheckout
 integer    :: jcount,jprint
!
!  define the constant zi sqrt(-1)
!
 zi=(-1.d-20,1.)
!
!  set up grid, extending from rin to rout using n gridpoints
!
 n=200
 rin=1.
 rout=90.
!
! define the sizes of non-Keplerian terms at Rin
!
 mode = 'blackhole'
 select case(mode)
 case('blackhole')
    call set_bh(spin=0.5585d0,rsch=0.5d0) ! Schwarzschild Radius: Rin = 2Rs
    call get_bh(rin,etazero,zetazero)
 case('binary')
    call set_binary(mass1=0.5d0,mass2=0.5d0,r1=0.25d0*rin,r2=0.25d0*rin)
    call get_binary(rin,etazero,zetazero)
 case default
    zetazero=0.
    etazero=0.      
 end select
 print*,' ETAZERO = ',etazero,' ZETAZERO = ',zetazero
!
!   define H/R at R=1
!
 honr=0.1030
!
! define the dissipation coefficient
!
 alpha=0.05

 write(6,"(1x, 'H/R, eta0, zeta0, alpha ', 4(es12.4))") honr, etazero, zetazero, alpha
 write(6,"(1x, 'have set zi = ', 2(es12.4))") zi
!
! define the position of the initial step in tiltangle
!
 rstep=20.
 wstep=2.

 write(6,"(1x, ' N, inner radius, outer radius', I6, 2(es12.4))") n, rin, rout
 write(6,"(1x, ' rstep   wstep  ', 2(es12.4))") rstep,wstep

 call makegrid  ! set up the radial grid
 call makedisc  ! set up the disc
 call setup     ! set up the initial conditions
!
! set up the run parameters
!
 jprint=10000000
 jcount=0
 nstep=0
 time=0.
 tstop=4000.
 tprint=2.*tstop
 tcheck=0.
 ctime=0.03

 toutfile=tstop/10.
 tcheckout=0.
 nfile=0
 write(6,"(1x, 'tstop toutfile ctime ', 3(es12.4))") tstop,toutfile,ctime
!
!   initial printout
!
 call prdisc
 call print
 call write_output_file
!
!   start the main evolution loop
!
!   calculate timestep - linear problem so do only once
!
 call tstep
 write(6,"(1x,'timestep dt = ', es12.4)") dt

 do while (time < tstop)
    nstep=nstep+1
    jcount=jcount+1
!
!  update the variables
!
    time=time+dt
    tcheck=tcheck+dt
    tcheckout=tcheckout+dt

    call update
!
!  print if necessary
!
    if(jcount.ge.jprint.or.tcheck.ge.tprint) then
       jcount=0
       tcheck=0.
       call print
    endif
!
!  write to output file if necessary
!
    if(tcheckout.ge.toutfile) then
       tcheckout=0.
       call write_output_file
    endif
 enddo

 call print
 call write_output_file

end program wave

!
! Subroutine to set up the grid points in r
!
subroutine makegrid
 use waveutils, only:n,r,rin,rout,n,dr,rsq,r12,r32
 implicit none
 integer :: i
 real*8  :: factor

 r(2)=rin
 r(2*n+2)=rout

 factor=(rout/rin)**(1./float(2*n))

 do i=3,2*n+1
    r(i)=rin*factor**(i-2)
 enddo

 do i=3,2*n+1
    dr(i)=r(i+1)-r(i-1)
 enddo

 dr(2)=dr(3)
 dr(2*n+2)=dr(2*n+1)

 do i=2,2*n+2
    rsq(i)=r(i)*r(i)
    r12(i)=sqrt(r(i))
    r32(i)=r(i)*r12(i)
 enddo

 return
end subroutine makegrid

!
! Subroutine to set up the disc parameters
! (surface density, temperature profile and orbital parameters)
!
subroutine makedisc
 use waveutils, only:n,r,r32,honr,csq,omega,eta,zeta,etazero,zetazero,rho
 implicit none
 integer :: i

 do i=2,2*n+2
    rho(i)=r(i)
!
!  note that rho = Sigma H^2
!     csq is suare of the sound speed
!     omega is the angular velocity of the disc
!
    csq(i)=((honr)**2)/r32(i)
    omega(i)=r(i)**(-1.5)

    eta(i)=etazero/r(i)
    zeta(i)=zetazero/r32(i)
 enddo

 return
end subroutine makedisc

!
!   this sets up the initial conditions in the disc
!
!   za1 and za2 correspond to A; zd1 and zd2 correspond to D
!   1 and 2 correspond to different levels of the Leapfrog.
!
subroutine setup
 use waveutils, only:n,r,pi,za1,za2,zd1,zd2,rstep,wstep,rsq
 implicit none
 integer :: i
 real*8  :: radius,tilt

 do i=1,n+1
    za1(i)=(0.,0.)
    za2(i)=(0.,0.)
    zd1(i)=(0.,0.)
    zd2(i)=(0.,0.)
 enddo

 do i=1,n
    radius=r(2*i+1)
    if(radius.lt.rstep-wstep) then
       tilt=0.
    elseif(radius.gt.rstep+wstep) then
       tilt=1.
    else
       tilt=0.5*(1.+sin(pi*(radius-rstep)/2./wstep))
    endif

    zd1(i)=cmplx(tilt/rsq(2*i+1),0.)
    zd2(i)=cmplx(tilt/rsq(2*i+1),0.)
 enddo

 return
end subroutine setup

!
!   Timestep constraints from viscosity parameters and orbital timescales
!
subroutine tstep
 use waveutils, only:n,dr,dt,csq,omega,eta,zi,ctime,zeta
 implicit none
 integer :: i,izone,itype
 real*8  :: tiny,dtry1,dtry2,dtry3
!
!  **** have not included dissipation term in this
!
!   if we add tidal forcing, then this should be included here too
!
 izone=0
 itype=0
 tiny=1.d-20
 dt=1.d20

 do i=2,n
    dtry1=dr(2*i)/sqrt(csq(2*i))
    dtry2=abs(zi)*omega(2*i)*eta(i)+tiny
    dtry2=0.1/(abs(dtry2))

    dtry3=abs(zi)*omega(2*i)*zeta(i)+tiny
    dtry3=0.1/(abs(dtry3))

    if(dtry1.lt.dt) then
       dt=dtry1
       izone=i
       itype=1
    endif

    if(dtry2.lt.dt) then
       dt=dtry2
       izone=i
       itype=2
    endif

    if(dtry3.lt.dt) then
       dt=dtry3
       izone=i
       itype=3
    endif

 enddo

 dt=ctime*dt

 if(itype.eq.1) write(6,"(1x, 'timestep due to wave speed at gridpt ',i6)") izone
 if(itype.eq.2) write(6,"(1x, 'timestep due eta at gridpt ',i6)") izone
 if(itype.eq.3) write(6,"(1x, 'timestep due zeta at gridpt ',i6)") izone

 return
end subroutine tstep

!
!--main bit: evolution of the equations through one timestep
!
subroutine update
 use waveutils
 implicit none
 integer :: i
!
! first we update zd1, which is at the half-gridpoints
!
 do i=1,n
    zd1(i)=zd1(i)-dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1) &
           *(rho(2*i+2)*r12(2*i+2)*za2(i+1)-rho(2*i)*r12(2*i)*za2(i))/dr(2*i+1) &
           +dt*zi*zeta(2*i+1)*omega(2*i+1)*zd2(i)
 enddo
!
! then update za1, which is at the full gridpoints:
! first keep za1=0 at the boundaries - for reflection off the boundaries
!
  za1(1)=(0.,0.)
  za1(n+1)=(0.,0.)
!
! then update the bulk of the grid
!
 do i=2,n
    za1(i)=za1(i) -dt*.5/rsq(2*i)*(rsq(2*i+1)*zd2(i)-rsq(2*i-1)*zd2(i-1))/dr(2*i) &
                  +dt*zi*eta(2*i)*omega(2*i)*za2(i) &
                  -dt*alpha*omega(2*i)*za1(i) 
 enddo
!
! then we update zd2, which is at the half-gridpoints
!
 do i=1,n
    zd2(i)=zd2(i) - dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1) &
                      *(rho(2*i+2)*r12(2*i+2)*za1(i+1)-rho(2*i)*r12(2*i)*za1(i))/dr(2*i+1) &
                  + dt*zi*zeta(2*i+1)*omega(2*i+1)*zd1(i)
 enddo
!
! then update za2, which is at the full gridpoints,
!
! first keep za2=0 at the boundaries
!
 za2(1)=(0.,0.)
 za2(n+1)=(0.,0.)
!
! then update the bulk of the grid
!
 do i=2,n
    za2(i)=za2(i) - dt*.5/rsq(2*i)*(rsq(2*i+1)*zd1(i)-rsq(2*i-1)*zd1(i-1))/dr(2*i) &
                  + dt*zi*eta(2*i)*omega(2*i)*za1(i) &
                  - dt*alpha*omega(2*i)*za2(i) 
 enddo

 return
end subroutine update

!
! Print info about the timestep
!
subroutine print
 use waveutils, only:n,r,zd1,za1,nstep,time
 implicit none
 integer :: i

 write(6,"('nstep= ',i8,' time= ',es12.4)") nstep, time

 do i=1,n+1
    write(6,"(1x, 'r,d,a  ', 5(es12.4))") r(2*i),zd1(i),za1(i)
 enddo

 return
end subroutine print

!
! Prints out the details of the grid and the disc
!
subroutine prdisc
 use waveutils, only:n,r,dr,rho,csq,omega,eta,zeta
 implicit none
 integer :: j
 
 write(6,"(1x, 'j,  r(j),  dr(j), rho(j), csq(j)')")
 do j=1,2*n+2
    write(6,"(1x, I4, 4(es12.4))") j,r(j),dr(j),rho(j),csq(j)
 enddo

 write(6,"(1x, 'j, r(j),  omega(j), eta(j), zeta(j)')")

 do j=1,2*n+2
    write(6,"(1x, I4, 4(es12.4))") j,r(j),omega(j),eta(j),zeta(j)
 enddo

 return
 end subroutine prdisc

!
! Prints output file
!
subroutine write_output_file
 use waveutils, only:n,rsq,zd1,zi,r,nstep,time,nfile
 implicit none
 integer    :: i
 complex*16 :: ztilt
 real*8     :: rtilt,xitilt,tilt,phase
 character(len=120) :: filename

 write(filename,"('wave_',i5.5)") nfile
 open(unit=24,file=filename,status='replace',form='formatted')
 write(6,"(a,i6,a,es12.4,a)") ' nstep ',nstep,' time = ',time,' writing '//trim(filename)

 write(24,*) time,nstep
 do i=1,n+1
    ztilt=zd1(i)*rsq(2*i+1)

    rtilt=ztilt
    xitilt=-zi*ztilt
    tilt=abs(ztilt)
    phase=atan2(xitilt,rtilt)

    write(24,"(1x,5F12.4)") r(2*i), ztilt, tilt, phase
 enddo
 close(unit=24)

 nfile = nfile + 1

 return
end subroutine write_output_file








