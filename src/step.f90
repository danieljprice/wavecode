module step
 implicit none
contains
!
!-- Timestep constraints from viscosity parameters and orbital timescales
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
 izone   = 0
 itype   = 0
 epsilon = tiny(epsilon)
 dt      = huge(dt)

 do i=2,n
    dtry1 = dr(2*i)/sqrt(csq(2*i))
    dtry2 = abs(zi)*omega(2*i)*eta(i)+epsilon
    dtry2 = 0.1/(abs(dtry2))

    dtry3 = abs(zi)*omega(2*i)*zeta(i)+epsilon
    dtry3 = 0.1/(abs(dtry3))

    if(dtry1.lt.dt) then
       dt    = dtry1
       izone = i
       itype = 1
    endif

    if(dtry2.lt.dt) then
       dt    = dtry2
       izone = i
       itype = 2
    endif

    if(dtry3.lt.dt) then
       dt    = dtry3
       izone = i
       itype = 3
    endif

 enddo

 dt = ctime*dt

 if(itype.eq.1) write(6,"(1x, 'timestep due to wave speed at gridpt ',i6)") izone
 if(itype.eq.2) write(6,"(1x, 'timestep due eta at gridpt ',i6)") izone
 if(itype.eq.3) write(6,"(1x, 'timestep due zeta at gridpt ',i6)") izone

 return
end subroutine tstep

!
!--main bit: evolution of the equations through one timestep
!
subroutine update
 use waveutils, only:n,dt,zi,zd1,zd2,za1,za2,omega,alpha,zeta,dr,eta,r12,rsq,rho,csq
 implicit none
 integer :: i

! first we update zd1, which is at the half-gridpoints
 do i=1,n
    zd1(i) = zd1(i)-dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1) &
              *(rho(2*i+2)*r12(2*i+2)*za2(i+1)-rho(2*i)*r12(2*i)*za2(i))/dr(2*i+1) &
              +dt*zi*zeta(2*i+1)*omega(2*i+1)*zd2(i)
 enddo

! then update za1, which is at the full gridpoints:
! first keep za1=0 at the boundaries - for reflection off the boundaries
 za1(1)   = (0.,0.)
 za1(n+1) = (0.,0.)

! then update the bulk of the grid
 do i=2,n
    za1(i) = za1(i) -dt*.5/rsq(2*i)*(rsq(2*i+1)*zd2(i)-rsq(2*i-1)*zd2(i-1))/dr(2*i) &
                     +dt*zi*eta(2*i)*omega(2*i)*za2(i) &
                     -dt*alpha(2*i)*omega(2*i)*za1(i)
 enddo

! then we update zd2, which is at the half-gridpoints
 do i=1,n
    zd2(i) = zd2(i) - dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1) &
                     *(rho(2*i+2)*r12(2*i+2)*za1(i+1)-rho(2*i)*r12(2*i)*za1(i))/dr(2*i+1) &
                     + dt*zi*zeta(2*i+1)*omega(2*i+1)*zd1(i)
 enddo
! then update za2, which is at the full gridpoints,

! first keep za2=0 at the boundaries
 za2(1)   = (0.,0.)
 za2(n+1) = (0.,0.)

! then update the bulk of the grid
 do i=2,n
    za2(i) = za2(i) - dt*.5/rsq(2*i)*(rsq(2*i+1)*zd1(i)-rsq(2*i-1)*zd1(i-1))/dr(2*i) &
                    + dt*zi*eta(2*i)*omega(2*i)*za1(i) &
                    - dt*alpha(2*i)*omega(2*i)*za2(i)
 enddo

 return
end subroutine update

end module step
