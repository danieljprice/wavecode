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
!----------------
!  MAIN CODE
!----------------
program wave
 use waveutils,       only:alphaSS,ctime,dt,etazero,zetazero,honr,mode,n,nfile,nstep
 use waveutils,       only:rin,rout,rstep,time,use_ext_sigma_profile,wstep,zi
 use waveutils,       only:use_ext_sigma_profile,p_index,q_index,use_pn,spin,theta
 use binary,          only:set_binary,get_binary
 use blackhole,       only:set_bh,get_bh
 use utils_setupfile, only:runtime_parameters

implicit none
 real     :: tcheck,tprint,tstop
 real     :: toutfile,tcheckout
 real     :: t1,t2,omegazero
 integer  :: jcount,jprint
 logical  :: iexist
 character(len=100) :: fname

!-- If file is passed as command argument, read sigma (& more) from it
 use_ext_sigma_profile = .false.
 if (command_argument_count()>0) then
    call get_command_argument(1,fname)
    inquire(file=fname,exist=iexist)
    if (iexist) then
       write(*,*) 'Using the sigma profile from: ',trim(fname)
    else
       write(*,*) "ERROR: The file '",trim(fname),"' does not exist"
       stop
    endif
    use_ext_sigma_profile = .true.
 endif

 zi=(0.0,1.) !  define the constant zi sqrt(-1)

!-- Default runtime parameters
 n       = 300         !  set up grid, extending from rin to rout using n gridpoints
 rin     = 10.
 rout    = 40.  !35.
 honr    = 0.05 !030   ! define H/R at R=1
 theta   = 3.          ! tile angle (degrees)
 mode    = 'blackhole' ! define the sizes of non-Keplerian terms at Rin

!--- Only used if mode=='blackhole'
 use_pn  = .false.
 spin    = 0.9

!--- Only used if .not.use_ext_sigma_profile
 p_index = 1.5
 q_index = 0.75
 alphaSS = 0.02        ! define the dissipation coefficient

 call runtime_parameters()

 select case(mode)
 case('blackhole')
    call set_bh(spin=spin,rsch=0.5) ! Schwarzschild Radius: Rin = 2Rs
    call get_bh(rin,etazero,zetazero,use_pn)
 case('binary')
    call set_binary(mass1=0.5,mass2=0.5,r1=0.25*rin,r2=0.25*rin)
    call get_binary(rin,etazero,zetazero,omegazero)
 case('binary-alpha0')
    call set_binary(mass1=0.5,mass2=0.5,r1=0.5*rin,r2=0.5*rin)
    call get_binary(rin,etazero,zetazero,omegazero)
    alphaSS = 0.2
    mode = 'binary'
 case default
    zetazero=0.
    etazero=0.
 end select
 print*,' ETAZERO = ',etazero,' ZETAZERO = ',zetazero,' OMEGAZERO = ',omegazero

 if (use_ext_sigma_profile) then
    write(6,"(1x, 'H/R, eta0, zeta0, alpha ', 3(es12.4))") honr, etazero, zetazero
 else
    write(6,"(1x, 'H/R, eta0, zeta0, alpha ', 4(es12.4))") honr, etazero, zetazero, alphaSS
 endif
 write(6,"(1x, 'have set zi = ', 2(es12.4))") zi

! define the position of the initial step in tiltangle
 rstep=20.
 wstep=2.

 write(6,"(1x, ' N, inner radius, outer radius', I6, 2(es12.4))") n, rin, rout
 write(6,"(1x, ' rstep   wstep  ', 2(es12.4))") rstep,wstep

 call makegrid  ! set up the radial grid
 call makedisc  ! set up the disc
 call setup     ! set up the initial conditions

! set up the run parameters
 jprint=10000000
 jcount=0
 nstep=0
 time=0.
 tstop=5000.!/8.+epsilon(0.)
 tprint=2.*tstop
 tcheck=0.
 ctime=0.03

 toutfile=tstop/40.
 tcheckout=0.
 nfile=0
 write(6,"(1x, 'tstop toutfile ctime ', 3(es12.4))") tstop,toutfile,ctime
!
!   initial printout
!
 !call prdisc
 !call print
 call write_output_file
!
!   start the main evolution loop
!
!   calculate timestep - linear problem so do only once
!
 call tstep
 write(6,"(1x,'timestep dt = ', es12.4)") dt
 call cpu_time(t1)

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
       !call print
    endif
!
!  write to output file if necessary
!
    if(tcheckout.ge.toutfile) then
       tcheckout=0.
       call write_output_file
       call cpu_time(t2)
       print "(a,f6.2)",' cpu time since last dump = ',t2 - t1
       t1 = t2
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
 real  :: factor,rinsav

 rinsav = rin
 rin = 1.

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

 rin = rinsav

 return
end subroutine makegrid

!
! Subroutine to set up the disc parameters
! (surface density, temperature profile and orbital parameters)
!
subroutine makedisc
 use waveutils, only:n,r,r32,honr,csq,omega,eta,zeta,etazero,zetazero,rho
 use waveutils, only:sigma,scale_height,use_ext_sigma_profile,alpha,alphaSS
 use waveutils, only:mode,p_index,q_index,use_pn
 use blackhole, only:get_bh
 use binary,    only:get_binary
 implicit none
 integer :: i,ierr,nlines,j
 integer, parameter :: isigma = 10
 real :: gradient, sigma_tolerance
 real, dimension(:), allocatable :: ext_sigma,ext_radius,ext_honh
 logical :: found_r
 real :: alphaAV
 character(len=100) :: fname

 sigma_tolerance = 1.e-14
 nlines = 0

 ! If a sigma profile is provided, read it in
  if (use_ext_sigma_profile) then
    call get_command_argument(1,fname)
    open(unit=isigma,file=fname,status='old',form='formatted',iostat=ierr)
    if (ierr/=0) STOP 'Could not open file!'

    write(*,'(a)') 'Please enter alphaAV used in the simulation: '
    read(*,*) alphaAV

    ! Work out how long the file is
    do while (ierr == 0)
      read(isigma,*,iostat=ierr)
      nlines = nlines + 1
    enddo
    nlines = nlines - 3 ! To take account of the header
    close(unit=isigma)

    ! Now save the profile
    allocate(ext_sigma(nlines),ext_radius(nlines),ext_honh(nlines))
    open(unit=isigma,file=fname,status='old',form='formatted',iostat=ierr)
    read(isigma,*)
    read(isigma,*)
    do i = 1,nlines
      read(isigma,*) ext_radius(i),ext_sigma(i),ext_honh(i)
    enddo
    close(unit=isigma)
    ext_radius = 0.25*ext_radius ! This is because of the scaling
  endif

 do i=2,2*n+2
!    rho(i) = Sigma*H**2 = (R**-p)*(R**(-q+3/2))^2
    if (use_ext_sigma_profile) then
      ! Linearly interpolate for each radial value what the sigma value should be
      found_r = .false.
      j = 1
      do while (.not.found_r)
        ! Find the correct radial bin
        if (ext_radius(j) <= r(i) .and. r(i) <= ext_radius(j+1)) then
          found_r = .true.
        else
          j = j + 1
          if (j > nlines) print*,'problem with provided sigma'
        endif
      enddo
      ! Interpolate between the nearest two points
      ! Split the equation into two lines for convenience
      gradient = (ext_sigma(j+1) - ext_sigma(j))/(ext_radius(j+1) - ext_radius(j))
      sigma(i) = ext_sigma(j) + (r(i) - ext_radius(j))*gradient
      if (abs(sigma(i)) < sigma_tolerance) sigma(i) = sigma_tolerance
      gradient = (ext_honh(j+1) - ext_honh(j))/(ext_radius(j+1) - ext_radius(j))
      alpha(i) = 1./10. * alphaAV * (ext_honh(j) + (r(i) - ext_radius(j))*gradient)
    else
      sigma(i)=r(i)**(-p_index) !*(1. - sqrt(1./r(i)))
      alpha(i)=alphaSS
    endif
    scale_height(i) = honr*r(i)**(-2*q_index + 3)
    rho(i) = sigma(i)*scale_height(i)
!
!  note that rho = Sigma H^2
!     csq is square of the sound speed
!     omega is the angular velocity of the disc
!
!    csq(i)=((honr)**2)/r32(i)
    csq(i)=((honr)**2)*r(i)**(-2.*q_index)
    omega(i)=r(i)**(-1.5)

    select case(trim(mode))
    case('blackhole')
       call get_bh(r(i),eta(i),zeta(i),use_pn)
    case('binary')
       call get_binary(r(i),eta(i),zeta(i),omega(i))
!       print*,' got BINARY',eta(i),zeta(i)
!       read*
    case default
       eta(i)=etazero/r(i)
       zeta(i)=zetazero/r32(i)
    end select
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
 use waveutils, only:n,r,pi,za1,za2,zd1,zd2,rstep,wstep,rsq,theta
 implicit none
 integer :: i
 real  :: radius,tilt

 do i=1,n+1
    za1(i)=(0.,0.)
    za2(i)=(0.,0.)
    zd1(i)=(0.,0.)
    zd2(i)=(0.,0.)
 enddo

 do i=1,n
    radius=r(2*i+1)
    ! if(radius.lt.rstep-wstep) then
    !    tilt=0.
    ! elseif(radius.gt.rstep+wstep) then
    !    tilt=1.
    ! else
    !    tilt=0.5*(1.+sin(pi*(radius-rstep)/2./wstep))
    ! endif
    tilt = sin(theta*pi/180.)

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
 real  :: epsilon,dtry1,dtry2,dtry3
!
!  **** have not included dissipation term in this
!
!   if we add tidal forcing, then this should be included here too
!
 izone=0
 itype=0
 epsilon=tiny(epsilon)
 dt=huge(dt)

 do i=2,n
    dtry1=dr(2*i)/sqrt(csq(2*i))
    dtry2=abs(zi)*omega(2*i)*eta(i)+epsilon
    dtry2=0.1/(abs(dtry2))

    dtry3=abs(zi)*omega(2*i)*zeta(i)+epsilon
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
                  -dt*alpha(2*i)*omega(2*i)*za1(i)
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
                  - dt*alpha(2*i)*omega(2*i)*za2(i)
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
 use waveutils, only:n,rsq,zd1,zi,r,nstep
 use waveutils, only:time,nfile,mode,sigma,alpha
 implicit none
 integer    :: i
 complex :: ztilt
 real     :: rtilt,xitilt,tilt,phase
 character(len=120) :: filename

 write(filename,"('angm',i5.5,'.exact')") nfile
 open(unit=24,file=filename,status='replace',form='formatted')
 write(6,"(a,i10,a,es12.4,a)") ' nstep ',nstep,' time = ',time,' writing '//trim(filename)

 if (trim(mode)=='blackhole') then
    print*,' time = ',time,' translating to ',8.*time
    write(24,*) 8.*time,nstep
 else
    write(24,*) time,nstep
 endif

 write(24,'(a)') '#  radius  sigma  ztilt  tilt  phase  alpha'

 do i=1,n+1
    ztilt=zd1(i)*rsq(2*i+1)

    rtilt=real(ztilt)
    xitilt=real(-zi*ztilt)
    tilt=abs(ztilt)
    phase=atan2(xitilt,rtilt)

    if (trim(mode)=='blackhole') then
       write(24,"(6(es18.10,1X))") 4.*r(2*i), sigma(2*i), real(ztilt), tilt, phase, alpha(2*i)
    else
       write(24,"(1x,5F12.4)") r(2*i), ztilt, tilt, phase
    endif
 enddo
 close(unit=24)

 nfile = nfile + 10

 return
end subroutine write_output_file
