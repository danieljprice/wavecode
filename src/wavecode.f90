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
 use output,          only:write_output_file,print_tstep
 use step,            only:update,tstep
 use setup,           only:makegrid,makedisc,do_setup
 use readwrite_infile,only:runtime_parameters

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
 call do_setup  ! set up the initial conditions

! set up the run parameters
 jprint = 10000000
 jcount = 0
 nstep  = 0
 time   = 0.
 tstop  = 5000.!/8.+epsilon(0.)
 tprint = 2.*tstop
 tcheck = 0.
 ctime  = 0.03

 toutfile  = tstop/40.
 tcheckout = 0.
 nfile     = 0
 write(6,"(1x, 'tstop toutfile ctime ', 3(es12.4))") tstop,toutfile,ctime
!
!   initial printout
!
 !call prdisc
 !call print_tstep
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
       !call print_tstep
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

 call print_tstep
 call write_output_file

end program wave
