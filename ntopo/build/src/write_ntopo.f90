module write_ntopo
USE nrtype
USE netcdf
implicit none
private
public::defineFile
public::write_iVec
public::write_dVec
! define dimension names
character(len=32),parameter :: sSeg_DimName='sSeg' ! dimension name for the stream segments
character(len=32),parameter :: sUps_DimName='sUps' ! dimension name for all upstream stream segments
character(len=32),parameter :: sAll_DimName='sAll' ! dimension name for total number of upstream reachs for all reaches 
character(len=32),parameter :: sHru_DimName='sHRU' ! dimension name for total number of upstream reachs for all reaches 
contains

 ! *********************************************************************
 ! new subroutine: define NetCDF file
 ! *********************************************************************
 subroutine defineFile(fname,           &  ! input: filename
                       nSeg,            &  ! input: number of stream segments
                       nTotal,          &  ! input: sum of number of upstream reaches for all reaches
                       nUps,            &  ! input: sum of number of immediate upstream reaches for all reaches
                       nHru,            &  ! input: sum of number of immediate upstream HRUs for all reaches
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 integer(i4b), intent(in)        :: nTotal       ! Sum of number of all upstream reaches for all reaches
 integer(i4b), intent(in)        :: nUps         ! Sum of number of immediate upstream reaches for all reaches
 integer(i4b), intent(in)        :: nHru         ! Sum of number of immediate upstream HRU for all reaches
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: idimID       ! dimension ID
 integer(i4b)                    :: iVar         ! variable index
 integer(i4b),parameter          :: nVars=30     ! number of variables
 integer(i4b),parameter          :: strLen=256   ! length of character string
 character(len=strLen)           :: cmessage     ! error message of downwind routine
 ! initialize error control
 ierr=0; message='defineFile/'

 ! define file
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create the stream segment dimension
 ierr = nf90_def_dim(ncid, trim(sSeg_DimName), nSeg, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for all upstream stream segments
 ierr = nf90_def_dim(ncid, trim(sUps_DimName), nUps, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for
 ierr = nf90_def_dim(ncid, trim(sAll_DimName), nTotal, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for
 ierr = nf90_def_dim(ncid, trim(sHru_DimName), nHru, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! define variables
 do iVar=1,nVars
  ! define variable
  select case(ivar)
   ! define network topology (integers)
   case( 1); call defvar('reachIndex',          'Reach Index (0,1,..nrch-1)',                       '-',     (/sSeg_DimName/), nf90_int,    ierr,cmessage)
   case( 2); call defvar('reachID',             'Reach ID',                                         '-',     (/sSeg_DimName/), nf90_int,    ierr,cmessage)
   ! define reach parameter (double precision)
   case( 3); call defvar('reachSlope',          'Slope of reach',                                   '-',     (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case( 4); call defvar('reachLength',         'Length of reach',                                  'm',     (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case( 5); call defvar('basinArea',           'Local basin area',                                 'm2',    (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case( 6); call defvar('upstreamArea',        'Area upstream of each reach',                      'm2',    (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case( 7); call defvar('totalArea',           'Basin area + Upstream area',                       'm2',    (/sSeg_DimName/), nf90_double, ierr,cmessage)
   ! define h2b parameter (double precision)
   case( 8); call defvar('hruIndex',           'Index of hru dimension',                            '-',     (/sHru_DimName/), nf90_int, ierr,cmessage)
   case( 9); call defvar('hru_id',             'Hru id',                                            '-',     (/sHru_DimName/), nf90_int, ierr,cmessage)
   case(10); call defvar('hru_lon',            'Longitude of hru centroid',                         'degree',(/sHru_DimName/), nf90_double, ierr,cmessage)
   case(11); call defvar('hru_lat',            'Latitude of hru centroid',                          'degree',(/sHru_DimName/), nf90_double, ierr,cmessage)
   case(12); call defvar('hru_elev',           'Average hru elevation',                             'm',     (/sHru_DimName/), nf90_double, ierr,cmessage)
   case(13); call defvar('hru_area',           'hru area',                                          'm2',    (/sHru_DimName/), nf90_double, ierr,cmessage)
   case(14); call defvar('hru_weight',         'Areal weight to total basin area',                  '-',     (/sHru_DimName/), nf90_double, ierr,cmessage)
   ! define network parameters (double precision)
   case(15); call defvar('reachLat1',          'Start latitude',                                    '-',     (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case(16); call defvar('reachLat2',          'End latitude',                                      '-',     (/sSeg_DimName/), nf90_double, ierr,cmessage)
   case(17); call defvar('reachLon1',          'Start longitude',                                   'degree',(/sSeg_DimName/), nf90_double, ierr,cmessage)
   case(18); call defvar('reachLon2',          'End longitude',                                     'degree',(/sSeg_DimName/), nf90_double, ierr,cmessage)
   case(19); call defvar('upReachTotalLength', 'Total upstream length',                             'm',     (/sAll_DimName/), nf90_double, ierr,cmessage)
   ! define network topology (integers)
   case(20); call defvar('downReachIndex',     'Immidiate Downstream reach index',                  '-',     (/sSeg_DimName/), nf90_int, ierr,cmessage)
   case(21); call defvar('downReachID',        'Immidiate Downstream reach ID',                     '-',     (/sSeg_DimName/), nf90_int, ierr,cmessage)
   case(22); call defvar('upReachIndex',       'Immidiate Upstream reach index',                    '-',     (/sUps_DimName/), nf90_int, ierr,cmessage)
   case(23); call defvar('upReachID',          'Immidiate Upstream reach ID',                       '-',     (/sUps_DimName/), nf90_int, ierr,cmessage)
   case(24); call defvar('reachList',          'List of all upstream reach indices',                '-',     (/sAll_DimName/), nf90_int, ierr,cmessage)
   case(25); call defvar('reachStart',         'start index for list of upstream reaches',          '-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)
   case(26); call defvar('reachCount',         'number of upstream reaches in each reach',          '-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)
   case(27); call defvar('upReachStart',       'start index for list of immediate upstream reaches','-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)
   case(28); call defvar('upReachCount',       'number of immediate upstream reaches in each reach','-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)
   case(29); call defvar('upHruStart',         'start index for list of upstream Hru',              '-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)
   case(30); call defvar('upHruCount',         'number of upstream Hru in each reach',              '-',     (/sSeg_DimName/), nf90_int,   ierr,cmessage)

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

end module write_ntopo
