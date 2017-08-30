module write_simoutput

USE nrtype
USE netcdf

implicit none

private

public::defineFile
public::defineStateFile
public::write_iVec
public::write_dVec

! define dimension names
character(len=32),parameter :: sSeg_DimName     ='sSeg'     ! dimension name for the stream segments
character(len=32),parameter :: sUps_Dimname     ='sUps'     ! dimension name for all upstream stream segments
character(len=32),parameter :: time_DimName     ='time'     ! dimension name for time
character(len=32),parameter :: sEns_DimName     ='sEns'     ! dimension name for ensembles (could be different runoff simululations)
character(len=32),parameter :: sWav_DimName     ='sWav'     ! dimension name for number of waves in stream segments 
character(len=32),parameter :: sTdh_DimName     ='sTdh'     ! dimension name for number of uh routed future overland flow in a stream segment 
character(len=32),parameter :: sTdh_irf_DimName ='sTdh_irf' ! dimension name for number of uh routed future channel flow in a stream segment
character(len=32),parameter :: sTB_DimName      ='sTB'      ! dimension name for number of time bound 
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
   case(11); call defvar('KWTroutedRunoff',   'KWT routed runoff in each reach',                  'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage)
   case(12); call defvar('UpBasRoutedRunoff', 'sum of upstream basin routed runoff in each reach','m3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage) ! added by NM
   case(13); call defvar('IRFroutedRunoff',   'IRF routed runoff in each reach',                  'm3/s',(/sSeg_DimName,time_Dimname/),nf90_double,ierr,cmessage) ! added by NM
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
 ! subroutine: define state NetCDF file
 ! *********************************************************************
 subroutine defineStateFile(fname,           &  ! input: filename
                            nEns,            &  ! input: number of ensemble
                            nSeg,            &  ! input: number of stream segments
                            ntdh,            &  ! input: number of future time steps for hill-slope UH routing 
                            ntdh_irf,        &  ! input: number of future time steps for irf river UH routing 
                            nMaxWave,        &  ! input: maixmum number of waves in a stream 
                            routOpt,         &  ! input: 0-> Both, 1-> IRF, 2 -> KWT, otherwise error
                            ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nEns         ! number of stream segments
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 integer(i4b), intent(in)        :: ntdh         ! number of future time steps for hill-slope UH routing 
 integer(i4b), intent(in)        :: ntdh_irf     ! number of future time steps for irf river UH routing 
 integer(i4b), intent(in)        :: nMaxWave     ! maximum number of waves in a stream segment 
 integer(i4b), intent(in)        :: routOpt      ! maximum number of waves in a stream segment 
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: idimID       ! dimension ID
 integer(i4b),parameter          :: strLen=256   ! length of character string
 character(len=strLen)           :: cmessage     ! error message of downwind routine

 ! initialize error control
 ierr=0; message='defineStateFile/'

 if (routOpt>2)then;ierr=10;message=trim(message)//'routOpt must be 0, 1, or 2';return;endif

 ! define file
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create the stream segment dimension
 ierr = nf90_def_dim(ncid, trim(sEns_DimName), nEns, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create the stream segment dimension
 ierr = nf90_def_dim(ncid, trim(sSeg_DimName), nSeg, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create time bound dimension
 ierr = nf90_def_dim(ncid, trim(sTB_DimName), 2, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for all upstream stream segments
 ierr = nf90_def_dim(ncid, trim(sTdh_DimName), ntdh, idimId)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! create dimension for future UH routed flow size
 if (routOpt==0 .or. routOpt==1)then
   ierr = nf90_def_dim(ncid, trim(sTdh_irf_DimName), ntdh_irf, idimId)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 endif
 ! create dimension for wave number in a reach 
 if (routOpt==0 .or. routOpt==2)then
   ierr = nf90_def_dim(ncid, trim(sWav_DimName), nMaxWave, idimId)
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 endif 

 ! define segment id
 call defvar('reachID',    'reach ID',                                    '-',   (/sSeg_DimName/),                              nf90_int,    ierr, cmessage)
 ! define time bound 
 call defvar('time_bound', 'time bound at last time step',                'sec', (/sTB_DimName/),                               nf90_double, ierr, cmessage)
 ! Hill-Slope 
 call defvar('QFUTURE',    'Hill slope convoluted flow',                  'm3/s',(/sSeg_DimName,sTdh_DimName,sEns_dimName/),    nf90_double, ierr, cmessage)
 ! IRF 
 if (routOpt==0 .or. routOpt==1)then
   call defvar('irfsize',    'number of future flow time step',           '-',   (/sSeg_DimName,sEns_dimName/),                 nf90_int,    ierr, cmessage)
   call defvar('QFUTURE_IRF','IRF convoluted river flow',                 'm3/s',(/sSeg_DimName,sTdh_irf_DimName,sEns_dimName/),nf90_double, ierr, cmessage)
 endif
 ! KWT 
 if (routOpt==0 .or. routOpt==2)then
   call defvar('wavesize',   'number of wave in a reach',                 '-',   (/sSeg_DimName,sEns_dimName/),                 nf90_int,    ierr, cmessage)
   call defvar('QF',         'flow of a wave ',                           'm3/s',(/sSeg_DimName,sWav_DimName,sEns_dimName/),    nf90_double, ierr, cmessage)
   call defvar('QM',         'modified flow of a wave',                   'm3/s',(/sSeg_DimName,sWav_DimName,sEns_dimName/),    nf90_double, ierr, cmessage)
   call defvar('TI',         'entry time of a wave in reach',             'sec', (/sSeg_DimName,sWav_DimName,sEns_dimName/),    nf90_double, ierr, cmessage)
   call defvar('TR',         'time when a wave is expected to exit reach','sec', (/sSeg_DimName,sWav_DimName,sEns_dimName/),    nf90_double, ierr, cmessage)
   call defvar('RF',         'routing flat (1 if wave has exited reach)', '-',   (/sSeg_DimName,sWav_DimName,sEns_dimName/),    nf90_int,    ierr, cmessage) 
 endif
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

  end subroutine

 end subroutine

 ! *********************************************************************
 ! subroutine: write an integer vector
 ! *********************************************************************
 subroutine write_iVec(fname,           &  ! input: filename
                       vname,           &  ! input: variable name
                       iVec,            &  ! input: variable data
                       iStart,          &  ! input: start index
                       iCount,          &  ! input: length of vector
                       ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname        ! variable name
 integer(i4b), intent(in)        :: iVec(:)      ! variable data
 integer(i4b), intent(in)        :: iStart(:)    ! start index
 integer(i4b), intent(in)        :: iCount(:)    ! length of vector
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
 ierr = nf90_put_var(ncid,iVarId,iVec,start=iStart,count=iCount)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close output file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine write_iVec

 ! *********************************************************************
 ! subroutine: write a double precision vector
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
