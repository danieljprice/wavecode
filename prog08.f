      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
c
c  this is prog08.f
c    .... with the precession put in by leapfrog 
c
c     this solves for the propagation of a warp wave through a Keplerian
c     disc. It uses the warp equations as given in Lubow, Ogilvie &
c     Pringle (2002) MNRAS 337, 706, Section 4, using the variables A and
c     D. For a Keplerian disc we set 'zeta' and 'eta' to be zero.
c
c     for eta and zeta zero, the precession via leapfrog is not
c     required. this means that much of the complexity can be removed.
c     we no longer need both za1 and za1, nor both zd1 and zd2.
c
c  we can also add dissipation in the form of a decay of 'A' the radial
c  velocity component, using 'alpha'
c
c 
c
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
      common/diss/alpha
      common/ics/rstep,wstep
c
c   define pi
c
      pi=4.*atan(1.)
c
c  define the constant zi sqrt(-1)
c
      zi=(-1.d-20,1.)
c
c   define the sizes of non-Keplerian terms at r1
c
      zetazero=0.
      etazero=0.
c
c   define H/R at R=1
c
      honr=0.1
c
c define the dissipation coefficient
c
      alpha=0.05
c
      write(6,122) honr, etazero, zetazero, alpha
 122  format(1x, 'H/R, eta0, zeta0, alpha ', 4(1pd12.4))
c
      write(6,121) zi
 121  format(1x, 'have set zi = ', 2(1pd12.4))
c
c  set up grid, extending from r1 to r2 using n gridpoints
c
      n=2000
      r1=1.
      r2=90.
c
c
c   define the position of the initial step in tiltangle
c
      rstep=20.
      wstep=2.
c
      write(6,123) n, r1, r2
 123  format(1x, ' N, inner radius, outer radius', I6, 2(1pd12.4))
c
      write(6,124)rstep,wstep
 124  format(1x, ' rstep   wstep  ', 2(1pd12.4))
c
      call makegrid
c
c
c   set up the disc
c
c
      call makedisc
c
c
c   set up inital conditions
c
c  we start with the disc tilted at a fixed angle
c
c
      call setup
c
c  set up the run parameters
c
      jprint=10000000
      jcount=0
c
      nstep=0
c
      time=0.
      tstop=4000.
      tprint=2.*tstop
      tcheck=0.
      ctime=0.03
c
      tmongo=tstop/10.
      tchmon=0.
      mline=1
c
c
      write(6,201) tstop,tmongo,ctime
 201  format(1x, 'tstop tmongo ctime ', 3(1pd12.4))
c
c   initial printout
c
      call prdisc
      call print
      call pmongo
c
c   start the main evolution loop
c
c   calculate timestep - linear problem so do only once
c
      call tstep
c
      write(6,110) dt
 110  format(1x,'timestep dt = ', 1pd12.4)
c
c
c
 10   continue
c
      nstep=nstep+1
      jcount=jcount+1
c
c
c  update the variables
c
      time=time+dt
      tcheck=tcheck+dt
      tchmon=tchmon+dt
c
      call update
c
c  print if necessary
c
      if(jcount.ge.jprint.or.tcheck.ge.tprint) then
         jcount=0
         tcheck=0.
         call print
      endif
c
c  output to mongo if necessary
c
      if(tchmon.ge.tmongo) then
         tchmon=0.
         call pmongo
      endif
c
      if(time.lt.tstop) then
         go to 10
      end if
c
c
      call print
      call pmongo
c
c
      stop
      end
c   
      subroutine makegrid
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
      r(2)=r1
      r(2*n+2)=r2
c
      factor=(r2/r1)**(1./float(2*n))
c
      do 17 i=3,2*n+1
         r(i)=r1*factor**(i-2)
 17   continue
c
      do 18 i=3,2*n+1
         dr(i)=r(i+1)-r(i-1)
 18   continue
c
      dr(2)=dr(3)
      dr(2*n+2)=dr(2*n+1)
c
      do 19 i=2,2*n+2
         rsq(i)=r(i)*r(i)
         r12(i)=sqrt(r(i))
         r32(i)=r(i)*r12(i)
 19   continue
c
      return
      end
c
c
c
      subroutine makedisc
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
c  
c
      do 21 i=2,2*n+2
         rho(i)=r(i)
c
c  note that rho = Sigma H^2
c     csq is suare of the sound speed
c     omega is the angular velocity of the disc
c
         csq(i)=((honr)**2)/r32(i)
         omega(i)=r(i)**(-1.5)
c
         eta(i)=etazero/r(i)
         zeta(i)=zetazero/r32(i)
c
 21   continue
c
c
      return
      end
c
      subroutine setup
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
      common/ics/rstep,wstep
c
c   define pi
c
      pi=4.*atan(1.)
c
c   this sets up the initial conditions in the disc
c
c   za1 and za2 correspond to A; zd1 and zd2 correspond to D
c   1 and 2 correspond to different levels of the Leapfrog.
c
      do 27 i=1,n+1
         za1(i)=(0.,0.)
         za2(i)=(0.,0.)
         zd1(i)=(0.,0.)
         zd2(i)=(0.,0.)
 27   continue
c
      do 28 i=1,n
         radius=r(2*i+1)
         if(radius.lt.rstep-wstep) then
            tilt=0.
         elseif(radius.gt.rstep+wstep) then
            tilt=1.
         else
            tilt=0.5*(1.+sin(pi*(radius-rstep)/2./wstep))
         endif
c
         zd1(i)=tilt/rsq(2*i+1)
         zd2(i)=tilt/rsq(2*i+1)
 28   continue
c
      return
      end
c
c
      subroutine tstep
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
c  **** have not included dissipation term in this
c
c   if we add tidal forcing, then this should be included here too
c
      izone=0
      itype=0
c
      tiny=1.d-20
c
      dt=1.d20
c
      do 31 i=2,n
c
         dtry1=dr(2*i)/sqrt(csq(2*i))
c
         dtry2=abs(zi)*omega(2*i)*eta(i)+tiny
         dtry2=0.1/(abs(dtry2))
c
         dtry3=abs(zi)*omega(2*i)*zeta(i)+tiny
         dtry3=0.1/(abs(dtry3))
c
c
         if(dtry1.lt.dt) then
            dt=dtry1
            izone=i
            itype=1
         endif
c
         if(dtry2.lt.dt) then
            dt=dtry2
            izone=i
            itype=2
         endif
c
         if(dtry3.lt.dt) then
            dt=dtry3
            izone=i
            itype=3
         endif
c
 31   continue
c
      dt=ctime*dt
c
c
      if(itype.eq.1) write(6,91) izone
 91   format(1x, 'timestep due to wave speed at gridpt ',i6)
      if(itype.eq.2) write(6,92) izone
 92   format(1x, 'timestep due eta at gridpt ',i6)
      if(itype.eq.3) write(6,93) izone
 93   format(1x, 'timestep due zeta at gridpt ',i6)
c
      return
      end
c
c
c
c
      subroutine update
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
      common/diss/alpha
c
c
c
c  first we update zd1, which is at the half-gridpoints
c
c
      do 1 i=1,n
         zd1(i)=zd1(i)-dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1)
     +     *(rho(2*i+2)*r12(2*i+2)*za2(i+1)-rho(2*i)*r12(2*i)*za2(i))
     +     /dr(2*i+1)
     +        +dt*zi*zeta(2*i+1)*omega(2*i+1)*zd2(i)
 1    continue
c
c
c   then update za1, which is at the full gridpoints,
c    
c
c  first keep za1=0 at the boundaries - for reflection off the boundaries
c
      za1(1)=(0.,0.)
      za1(n+1)=(0.,0.)
c
c   then update the bulk of the grid
c
      do 2 i=2,n
         za1(i)=za1(i)
     +     -dt*.5/rsq(2*i)*(rsq(2*i+1)*zd2(i)-rsq(2*i-1)*zd2(i-1))
     +               /dr(2*i)
     +           +dt*zi*eta(2*i)*omega(2*i)*za2(i)
     +           -alpha*omega(2*i)*dt*za1(i) 
2     continue
c
c  then we update zd2, which is at the half-gridpoints
c
c
      do 11 i=1,n
         zd2(i)=zd2(i)-dt*csq(2*i+1)/2./r12(2*i+1)/rho(2*i+1)
     +     *(rho(2*i+2)*r12(2*i+2)*za1(i+1)-rho(2*i)*r12(2*i)*za1(i))
     +     /dr(2*i+1)
     +        +dt*zi*zeta(2*i+1)*omega(2*i+1)*zd1(i)
 11   continue
c
c
c   then update za2, which is at the full gridpoints,
c
c  first keep za2=0 at the boundaries
c
      za2(1)=(0.,0.)
      za2(n+1)=(0.,0.)
c
c   then update the bulk of the grid
c
      do 21 i=2,n
         za2(i)=za2(i)
     +     -dt*.5/rsq(2*i)*(rsq(2*i+1)*zd1(i)-rsq(2*i-1)*zd1(i-1))
     +               /dr(2*i)
     +           +dt*zi*eta(2*i)*omega(2*i)*za1(i)
     +           -alpha*omega(2*i)*dt*za2(i) 
 21   continue
c
c
      return
      end
c
c
      subroutine print
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
c
      write(6,101) nstep, time
c
      do 50 i=1,n+1
         write(6,102) r(2*i),zd1(i),za1(i)
 50   continue
c
c
 101  format(1x, 'nstep= ',i8,' time= ',1pd12.4)
 102  format(1x, 'r,d,a  ', 5(1pd12.4))
c
c
      return
      end
c
c
      subroutine prdisc
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
c  prints out the details of the grid and the disc
c
      write(6,140)
c
      do 41 j=1,2*n+2
         write(6,141) j,r(j),dr(j),rho(j),csq(j)
 41   continue
c
 140  format(1x, 'j,  r(j),  dr(j), rho(j), csq(j)')
 141  format(1x, I4, 4(1pd12.4))
c
      write(6,142)
c
      do 43 j=1,2*n+2
         write(6,143) j,r(j),omega(j),eta(j),zeta(j)
 43   continue
c
 142  format(1x, 'j, r(j),  omega(j), eta(j), zeta(j)')
 143  format(1x, I4, 4(1pd12.4))
c
      return
      end
c
c
      subroutine pmongo
      implicit real*8(a-h,o-z)
      complex*16 za1,za2,zd1,zd2 ,zi
      complex*16 ztilt
      common/zvbles/za1(2002),za2(2002),zd1(2002),zd2(2002),zi
      common/grid/r(4004),dr(4004),rsq(4004),
     +   r12(4004),r32(4004),r1,r2,n
      common/disc/omega(4004),rho(4004),csq(4004)
      common/tempus/time,dt,ctime,nstep,mline
      common/precess/eta(4004),zeta(4004)
      common/consts/honr,etazero,zetazero
c
c  prints where we are in format for mongo
c
      mline1=mline
      mline2=mline1+n
c
      write(6,150) mline1,mline2,nstep,time
 150  format(1x, 'lines', i6,' to',i6,' are nstep',i6,' time=', 1pd12.4)
c
      mline=mline2+1
c
      do 51 i=1,n+1
         ztilt=zd1(i)*rsq(2*i+1)
c
         rtilt=ztilt
         xitilt=-zi*ztilt
         tilt=abs(ztilt)
         phase=atan2(xitilt,rtilt)
c
         write(24,151) r(2*i), ztilt, tilt, phase
 51   continue
c
 151  format(1x, 5F12.4)
c
      return
      end








