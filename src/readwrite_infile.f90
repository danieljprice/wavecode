module readwrite_infile
 use waveutils, only:n,rin,rout,honr,alphaSS,mode,p_index,q_index
 use waveutils, only:use_ext_sigma_profile,spin,theta,rstep,wstep,tstop,nout
 use blackhole, only:ifreq
 use setup,     only:iwarp

 implicit none

 private

 public :: runtime_parameters

contains
!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# Setup file for wavecode:'
 call write_inopt(tstop  ,'tstop' ,'time to integrate to'                   ,iunit)
 call write_inopt(nout   ,'nout'  ,'number of output files'                 ,iunit)
 call write_inopt(n      ,'n'     ,'number of grid points'                  ,iunit)
 call write_inopt(honr   ,'honr'  ,'scale height H/R of disc (at Rin)'      ,iunit)
 call write_inopt(mode   ,'mode'  ,'blackhole, binary, or binar-alpha0'     ,iunit)
 call write_inopt(iwarp  ,'iwarp' ,'1=constant tilt angle, 2=tilt with step',iunit)

 if (iwarp==1) then
    write(iunit,"(/a)") '#------ Only used if iwarp == 1 -----------------------------------------------------'
    call write_inopt(theta,'theta','inclination of disc (degrees)'          ,iunit)
 endif

 if (iwarp==2) then
    write(iunit,"(/a)") '#------ Only used if iwarp == 2 -----------------------------------------------------'
    call write_inopt(rstep,'rstep','radius where step occurs (scaled by rin)',iunit)
    call write_inopt(wstep,'wstep','width of step (scaled by rin)',iunit)
 endif

 if (trim(mode)=='blackhole') then
    write(iunit,"(/a)") '#------ Only used if mode == blackhole ----------------------------------------------'
    call write_inopt(ifreq ,'ifreq','BH frequencies (1=GR, 2=GR-approx, 3=PN)',iunit)
    call write_inopt(spin  ,'spin' ,'spin parameter of black hole |a|<1'      ,iunit)
 endif

 if (.not.use_ext_sigma_profile) then
    write(iunit,"(/a)") '#------ Only used if not reading sigma (& more) from an external file ---------------'
    call write_inopt(rin    ,'rin'    ,'inner edge'                        ,iunit)
    call write_inopt(rout   ,'rout'   ,'outer edge'                        ,iunit)
    call write_inopt(alphaSS,'alphaSS','dissipation parameter'                     ,iunit)
    call write_inopt(p_index,'p_index','power law index of surface density profile',iunit)
    call write_inopt(q_index,'q_index','power law index of sound speed profile'    ,iunit)
 endif

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(tstop   ,'tstop',db,min=0.,      errcount=nerr)
 call read_inopt(nout    ,'nout' ,db,min=1 ,      errcount=nerr)
 call read_inopt(n       ,'n'    ,db,min=0 ,      errcount=nerr)
 call read_inopt(honr    ,'honr' ,db,min=0.,      errcount=nerr)
 call read_inopt(mode    ,'mode' ,db,             errcount=nerr)
 call read_inopt(iwarp   ,'iwarp',db,min=1 ,max=2,errcount=nerr)

 if (.not.use_ext_sigma_profile) then
    call read_inopt(rin    ,'rin'    ,db,min=0.,errcount=nerr)
    call read_inopt(rout   ,'rout'   ,db,min=0.,errcount=nerr)
    call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
    call read_inopt(p_index,'p_index',db,       errcount=nerr)
    call read_inopt(q_index,'q_index',db,       errcount=nerr)
 endif

 if (iwarp==1) then
    call read_inopt(theta,'theta',db,min=0.,max=90.,errcount=nerr)
 endif

 if (iwarp==2) then
    call read_inopt(rstep ,'rstep',db,min=0.,errcount=nerr)
    call read_inopt(wstep ,'wstep',db,min=0.,errcount=nerr)
 endif

 if (trim(mode)=='blackhole') then
    call read_inopt(ifreq ,'ifreq',db,min=1  ,max=3 ,errcount=nerr)
    call read_inopt(spin  ,'spin' ,db,min=-1.,max=1.,errcount=nerr)
 endif

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-- Read runtime parameters from setup file
subroutine runtime_parameters()
 character(*), parameter :: filename = 'setup.in'
 logical :: iexist
 integer :: ierr
 integer :: imode

 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then

    print*,'Please pick a mode:'
    print*,' 1 = blackhole'
    print*,' 2 = binary'
    print*,' 3 = binary-alpha0'
    read*,imode

    select case(imode)
    case(1)
       mode = 'blackhole'
    case(2)
       mode = 'binary'
    case(3)
       mode = 'binary-alpha0'
       alphaSS = 0.2
    case default
       stop 'ERROR: Bad choice of imode'
    end select

    call write_setupfile(filename)
    print*,' Edit '//trim(filename)//' and rerun wavecode'
    stop
 endif

end subroutine runtime_parameters

end module readwrite_infile
