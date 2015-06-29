module write_simoutput
USE nrtype
USE netcdf
implicit none
private
public::defineFile
public::write_iVec
public::write_dVec
! define dimension names
character(len=32),parameter :: sSeg_DimName='sSeg' ! dimension name for the stream segments
character(len=32),parameter :: sUps_Dimname='sUps' ! dimension name for all upstream stream segments
character(len=32),parameter :: time_DimName='time' ! dimension name for time
contains

 ! *********************************************************************
 ! new subroutine: define NetCDF file
 ! *********************************************************************
 subroutine defineFile(fname,           &  ! input: filename
                       nSeg,            &  ! input: number of stream segments
                       nTotal,          &  ! input: total number of upstream reaches for all reaches
                       units_time,      &  ! input: time units
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 integer(i4b), intent(in)        :: nTotal       ! total number of upstream reaches for all reaches
 character(*), intent(in)        :: units_time   ! time units
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: idimID       ! dimension ID
 integer(i4b)                    :: iVar         ! variable index
 integer(i4b),parameter          :: nVars=13     ! number of variables
 integer(i4b),parameter          :: strLen=256   ! length of character string
 character(len=strLen)           :: cmessage     ! error message of downwind routine
 ! initialize error control
 ierr=0; message='defineFile/'

 ! define file
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create time dimension (unlimited)
 ierr = nf90_def_dim(ncid, trim(time_DimName), nf90_unlimited, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create the stream segment dimension
 ierr = nf90_def_dim(ncid, trim(sSeg_DimName), nSeg, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for all upstream stream segments
 ierr = nf90_def_dim(ncid, trim(sUps_DimName), nTotal, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! define coordinate variable for time
 call defvar(trim(time_DimName),trim(time_DimName),trim(units_time),(/time_DimName/),nf90_double,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! define variables
 do iVar=1,nVars
  ! define variable
  select case(ivar)
   ! define network topology (integers)
   case( 1); call defvar('reachID',           'reach ID',                                '-',   (/sSeg_DimName/),             nf90_int,   ierr,cmessage)
   case( 2); call defvar('reachOrder',        'processing order',                        '-',   (/sSeg_DimName/),             nf90_int,   ierr,cmessage)
   case( 3); call defvar('reachList',         'list of upstream reaches',                '-',   (/sUps_DimName/),             nf90_int,   ierr,cmessage)
   case( 4); call defvar('listStart',         'start index for list of upstream reaches','-',   (/sSeg_DimName/),             nf90_int,   ierr,cmessage)
   case( 5); call defvar('listCount',         'number of upstream reaches in each reach','-',   (/sSeg_DimName/),             nf90_int,   ierr,cmessage)
   ! define network parameters (double precision)
   case( 6); call defvar('basinArea',         'local basin area',                        'm2',  (/sSeg_DimName/),             nf90_double,ierr,cmessage)
   case( 7); call defvar('upstreamArea',      'area upstream of each reach',             'm2',  (/sSeg_DimName/),             nf90_double,ierr,cmessage)
   ! define runoff variables
   case( 8); call defvar('instBasinRunoff',   'instantaneous basin runoff in each reach',         'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage)
   case( 9); call defvar('dlayBasinRunoff',   'delayed basin runoff in each reach',               'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage)
   case(10); call defvar('sumUpstreamRunoff', 'sum of upstream runoff in each reach',             'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage)
   case(11); call defvar('routedRunoff',      'routed runoff in each reach',                      'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage)
   case(12); call defvar('UpBasRoutedRunoff', 'sum of upstream basin routed runoff in each reach','m3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage) ! added by NM
   case(13); call defvar('VICroutedRunoff',   'VIC routed runoff in each reach',                  'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage) ! added by NM
   case default; ierr=20; message=trim(message)//'unable to identify variable index'; return
  end select
  ! check errors
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do

 ! end definitions
 ierr = nf90_enddef(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif


 contains
  subroutine defvar(vname,vdesc,vunit,dimNames,ivtype,ierr,message)
  ! input
  character(*), intent(in)   :: vname       ! variable name
  character(*), intent(in)   :: vdesc       ! variable description
  character(*), intent(in)   :: vunit       ! variable units
  character(*), intent(in)   :: dimNames(:) ! variable dimension names
  integer(i4b), intent(in)   :: ivtype      ! variable type
  ! output
  integer(i4b), intent(out)  :: ierr        ! error code
  character(*), intent(out)  :: message     ! error message
  ! local
  integer(i4b)               :: id          ! loop through dimensions
  integer(i4b)               :: dimIDs(size(dimNames))  ! vector of dimension IDs
  integer(i4b)               :: iVarId      ! variable ID

  ! define dimension IDs
  do id=1,size(dimNames)
   ierr=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end do

  ! define variable
  ierr = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! add variable description
  ierr = nf90_put_att(ncid,iVarId,'long_name',trim(vdesc))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  ! add variable units
  ierr = nf90_put_att(ncid,iVarId,'units',trim(vunit))
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  end subroutine defvar

 end subroutine defineFile

 ! *********************************************************************
 ! new subroutine: write an integer vector
 ! *********************************************************************
 subroutine write_iVec(fname,           &  ! input: filename
                       vname,           &  ! input: variable name
                       iVec,            &  ! input: variable data
                       iStart,          &  ! input: start index
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname        ! variable name
 integer(i4b), intent(in)        :: iVec(:)      ! variable data
 integer(i4b), intent(in)        :: iStart       ! start index
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarId       ! NetCDF variable ID
 ! initialize error control
 ierr=0; message='write_iVec/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_write,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! write data
 ierr = nf90_put_var(ncid,iVarId,iVec,start=(/iStart/),count=(/size(iVec)/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine write_iVec



 ! *********************************************************************
 ! new subroutine: write a double precision vector
 ! *********************************************************************
 subroutine write_dVec(fname,           &  ! input: filename
                       vname,           &  ! input: variable name
                       dVec,            &  ! input: variable data
                       iStart,          &  ! input: start index
                       iCount,          &  ! input: length of vector
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname        ! variable name
 real(dp), intent(in)            :: dVec(:)      ! variable data
 integer(i4b), intent(in)        :: iStart(:)    ! start indices
 integer(i4b), intent(in)        :: iCount(:)    ! length of vector
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarId       ! NetCDF variable ID
 ! initialize error control
 ierr=0; message='write_dVec/'

 ! open NetCDF file
 ierr = nf90_open(trim(fname),nf90_write,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get variable ID
 ierr = nf90_inq_varid(ncid,trim(vname),iVarId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! write data
 ierr = nf90_put_var(ncid,iVarId,dVec,start=iStart,count=iCount)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine write_dVec



end module write_simoutput
