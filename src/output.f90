module output
 implicit none
contains

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

!
! Print info about the timestep
!
subroutine print_tstep
 use waveutils, only:n,r,zd1,za1,nstep,time
 implicit none
 integer :: i

 write(6,"('nstep= ',i8,' time= ',es12.4)") nstep, time

 do i=1,n+1
    write(6,"(1x, 'r,d,a  ', 5(es12.4))") r(2*i),zd1(i),za1(i)
 enddo

 return
end subroutine print_tstep

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

end module output
