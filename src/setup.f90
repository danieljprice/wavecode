module setup
 implicit none

!--- Variables when reading sigma from an external file
 real, dimension(:), allocatable :: ext_sigma,ext_radius,ext_honh
 real    :: alphaAV
 integer :: nlines

contains
!
! Subroutine to set up the grid points in r
!
subroutine makegrid
 use waveutils, only:n,r,rin,rout,n,dr,rsq,r12,r32
 integer :: i
 real    :: factor,rinsav,routsav

 rinsav   = rin
 routsav  = rout
 rout     = rout/rin
 rin      = 1.

 r(2)     = rin
 r(2*n+2) = rout
 factor   = (rout/rin)**(1./float(2*n))

 do i=3,2*n+1
    r(i) = rin*factor**(i-2)
 enddo

 do i=3,2*n+1
    dr(i) = r(i+1)-r(i-1)
 enddo

 dr(2)     = dr(3)
 dr(2*n+2) = dr(2*n+1)

 do i=2,2*n+2
    rsq(i) = r(i)*r(i)
    r12(i) = sqrt(r(i))
    r32(i) = r(i)*r12(i)
 enddo

 rin = rinsav
 rout = routsav

 return
end subroutine makegrid

 !
 ! Subroutine to set up the disc parameters
 ! (surface density, temperature profile and orbital parameters)
 !
subroutine makedisc
 use waveutils, only:n,r,r32,honr,csq,omega,eta,zeta,etazero,zetazero,rho
 use waveutils, only:sigma,scale_height,use_ext_sigma_profile,alpha,alphaSS
 use waveutils, only:mode,p_index,q_index,rin
 use blackhole, only:get_bh
 use binary,    only:get_binary
 integer :: i,j
 real    :: gradient,sigma_tolerance
 logical :: found_r

 sigma_tolerance = 1.e-14
 if (use_ext_sigma_profile) then ! rescale radius by rin (which should be read from file)
    ext_radius = ext_radius/rin
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
             if (j > nlines) then
                 print*,'ERROR: problem with provided sigma file'
                 print*,'Could not find radial bin in file for r = ',r(i)
                 STOP
             endif
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
       sigma(i) = r(i)**(-p_index) !*(1. - sqrt(1./r(i)))
       alpha(i) = alphaSS
    endif
    scale_height(i) = honr*r(i)**(-2*q_index + 3) ! ?????
    rho(i) = sigma(i)*scale_height(i)             ! Does this get done right if sigma is read from file?
    !
    !  note that rho = Sigma H^2
    !     csq is square of the sound speed
    !     omega is the angular velocity of the disc
    !
    !    csq(i)=((honr)**2)/r32(i)
    csq(i)   = ((honr)**2)*r(i)**(-2.*q_index) !?????
    omega(i) = r(i)**(-1.5)                    !?????

    select case(trim(mode))
    case('blackhole')
       call get_bh(r(i),eta(i),zeta(i))
    case('binary')
       call get_binary(r(i),eta(i),zeta(i),omega(i))
       !       print*,' got BINARY',eta(i),zeta(i)
       !       read*
    case default
       eta(i)  = etazero/r(i)
       zeta(i) = zetazero/r32(i)
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
subroutine do_setup
 use waveutils, only:n,r,pi,za1,za2,zd1,zd2,rstep,wstep,rsq,theta
 integer :: i
 real  :: radius,tilt

 do i=1,n+1
    za1(i) = (0.,0.)
    za2(i) = (0.,0.)
    zd1(i) = (0.,0.)
    zd2(i) = (0.,0.)
 enddo

 do i=1,n
    radius=r(2*i+1)
    ! if(radius.lt.rstep-wstep) then
    !    tilt = 0.
    ! elseif(radius.gt.rstep+wstep) then
    !    tilt = 1.
    ! else
    !    tilt = 0.5*(1.+sin(pi*(radius-rstep)/2./wstep))
    ! endif
    tilt   = sin(theta*pi/180.)

    zd1(i) = cmplx(tilt/rsq(2*i+1),0.)
    zd2(i) = cmplx(tilt/rsq(2*i+1),0.)
 enddo

 return
end subroutine do_setup

subroutine read_external_sigma
 integer, parameter :: isigma = 10
 character(len=100) :: fname
 integer :: i,ierr

 nlines = 0

! If a sigma profile is provided, read it in
 call get_command_argument(1,fname)
 open(unit=isigma,file=fname,status='old',form='formatted',iostat=ierr)
 if (ierr/=0) STOP 'Could not open file!'

 write(*,'(a)') ' Please enter alphaAV used in the simulation: '
 read(*,*) alphaAV
 print*

 ! Work out how long the file is
 do while (ierr == 0)
    read(isigma,*,iostat=ierr)
    nlines = nlines + 1
 enddo
 nlines = nlines - 3 ! To take account of the header (2 lines) and because we added one to counter at end of file
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

end subroutine read_external_sigma

end module setup
