module utils_setupfile
 use waveutils, only:n,rin,rout,honr,alphaSS,mode

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
 write(iunit,"(a)") '# setup file for wavecode'
 call write_inopt(n      ,'n'      ,'number of grid points'             ,iunit)
 call write_inopt(rin    ,'rin'    ,'inner edge'                        ,iunit)
 call write_inopt(rout   ,'rout'   ,'outer edge'                        ,iunit)
 call write_inopt(honr   ,'honr'   ,'scale height H/R of disc (at R=1)' ,iunit)
 call write_inopt(alphaSS,'alphaSS','dissipation parameter'             ,iunit)
 call write_inopt(mode   ,'mode'   ,'blackhole, binary, or binar-alpha0',iunit)
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
 call read_inopt(n      ,'n'      ,db,min=0 ,errcount=nerr)
 call read_inopt(rin    ,'rin'    ,db,min=0.,errcount=nerr)
 call read_inopt(rout   ,'rout'   ,db,min=0.,errcount=nerr)
 call read_inopt(honr   ,'honr'   ,db,min=0.,errcount=nerr)
 call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 call read_inopt(mode   ,'mode'   ,db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-- Read runtime parameters from setup file
subroutine runtime_parameters()
 character(*), parameter :: filename = 'setup.file'
 logical :: iexist
 integer :: ierr

 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    call write_setupfile(filename)
    print*,' Edit '//trim(filename)//' and rerun wavecode'
    stop
 endif

end subroutine runtime_parameters

end module utils_setupfile
