MODULE write_simoutput

! Moudle wide external modules
USE nrtype
USE netcdf
USE public_var
implicit none

private

public::defineFile

CONTAINS

 ! *********************************************************************
 ! new subroutine: define routing output NetCDF file
 ! *********************************************************************
 SUBROUTINE defineFile(fname,           &  ! input: filename
                       nEns,            &  ! input: number of ensembles
                       nHRU,            &  ! input: number of HRUs
                       nSeg,            &  ! input: number of stream segments
                       units_time,      &  ! input: time units
                       calendar,        &  ! input: calendar
                       ierr, message)      ! output: error control
 !Dependent modules
 USE globalData, ONLY: meta_qDims
 USE var_lookup, ONLY: ixQdims, nQdims

 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 integer(i4b), intent(in)        :: nEns         ! number of ensembles
 integer(i4b), intent(in)        :: nHRU         ! number of HRUs
 integer(i4b), intent(in)        :: nSeg         ! number of stream segments
 character(*), intent(in)        :: units_time   ! time units
 character(*), intent(in)        :: calendar     ! calendar
 ! output variables
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: jDim, iVar   ! dimension, and variable index
 integer(i4b),parameter          :: nVars=8      ! number of variables
 character(len=strLen)           :: cmessage     ! error message of downwind routine

 ! initialize error control
 ierr=0; message='defineFile/'

 associate (dim_seg  => meta_qDims(ixQdims%seg)%dimName,    &
            dim_hru  => meta_qDims(ixQdims%hru)%dimName,    &
            dim_ens  => meta_qDims(ixQdims%ens)%dimName,    &
            dim_time => meta_qDims(ixQdims%time)%dimName)

! populate q dimension meta (not sure if this should be done here...)
 meta_qDims(ixQdims%seg)%dimLength = nSeg
 meta_qDims(ixQdims%hru)%dimLength = nHRU
 meta_qDims(ixQdims%ens)%dimLength = nEns

 ! --------------------
 ! define file
 ! --------------------
 ierr = nf90_create(trim(fname),nf90_classic_model,ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 do jDim =1,nQdims
   if (jDim ==ixQdims%time) then ! time dimension (unlimited)
    ierr = nf90_def_dim(ncid, trim(meta_qDims(jDim)%dimName), nf90_unlimited, meta_qDims(jDim)%dimId)
   else
    ierr = nf90_def_dim(ncid, trim(meta_qDims(jDim)%dimName), meta_qDims(jDim)%dimLength ,meta_qDims(jDim)%dimId)
   endif
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do

 ! define coordinate variable for time
 call defvar(ncid, trim(dim_time), (/dim_time/), nf90_double, ierr, cmessage, vdesc=trim(dim_time), vunit=trim(units_time), vcal=calendar)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! define variables
 do iVar=1,nVars
  ! define variable
  select case(iVar)
   ! define network topology (integers)
   case( 1); call defvar(ncid, 'basinID',           (/dim_hru/),          nf90_int,   ierr,cmessage, vdesc='basin ID',                            vunit='-'   )
   case( 2); call defvar(ncid, 'reachID',           (/dim_seg/),          nf90_int,   ierr,cmessage, vdesc='reach ID',                            vunit='-'   )
   ! define runoff variables (double precision)
   case( 3); call defvar(ncid, 'basRunoff',         (/dim_hru,dim_time/), nf90_float, ierr,cmessage, vdesc='basin runoff',                        vunit='m/s' )
   case( 4); call defvar(ncid, 'instRunoff',        (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='instantaneous runoff in each reach',  vunit='m3/s')
   case( 5); call defvar(ncid, 'dlayRunoff',        (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='delayed runoff in each reach',        vunit='m3/s')
   case( 6); call defvar(ncid, 'sumUpstreamRunoff', (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='sum of upstream runoff in each reach',vunit='m3/s')
   case( 7); call defvar(ncid, 'KWTroutedRunoff',   (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='KWT routed runoff in each reach',     vunit='m3/s')
   case( 8); call defvar(ncid, 'IRFroutedRunoff',   (/dim_seg,dim_time/), nf90_float, ierr,cmessage, vdesc='IRF routed runoff in each reach',     vunit='m3/s')
   case default; ierr=20; message=trim(message)//'unable to identify variable index'; return
  end select
  ! check errors
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 end do

 end associate

 ! end definitions
 ierr = nf90_enddef(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 END SUBROUTINE defineFile

 ! *********************************************************************
 ! private subroutine: define variable attributes NetCDF file
 ! *********************************************************************
 SUBROUTINE defvar(ncid, vname, dimNames, ivtype, ierr, message, vdesc, vunit, vcal)
  ! input
  integer(i4b), intent(in)             :: ncid                   ! Input: netcdf fine ID
  character(*), intent(in)             :: vname                  ! Input: variable name
  character(*), intent(in)             :: dimNames(:)            ! Input: variable dimension names
  integer(i4b), intent(in)             :: ivtype                 ! Input: variable type
  character(*), intent(in), optional   :: vdesc                  ! Input: variable description
  character(*), intent(in), optional   :: vunit                  ! Input: variable units
  character(*), intent(in), optional   :: vcal                   ! Input: calendar (if time variable)
  ! output
  integer(i4b), intent(out)            :: ierr                   ! error code
  character(*), intent(out)            :: message                ! error message
  ! local
  character(len=strLen)                :: calendar_str           ! calendar string
  character(len=strLen)                :: unit_str               ! unit string
  character(len=strLen)                :: desc_str               ! long_name string
  integer(i4b)                         :: id                     ! loop through dimensions
  integer(i4b)                         :: dimIDs(size(dimNames)) ! vector of dimension IDs
  integer(i4b)                         :: iVarId                 ! variable ID

  ! define dimension IDs
  do id=1,size(dimNames)
   ierr=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id))
   if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end do

  ! define variable
  ierr = nf90_def_var(ncid,trim(vname),ivtype,dimIds,iVarId)
  if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

  if (present(vdesc)) then ! add long_name
    desc_str = trim(vdesc)
    ierr = nf90_put_att(ncid,iVarId,'long_name',trim(desc_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

  if (present(vunit)) then ! add variable unit
    unit_str = trim(vunit)
    ierr = nf90_put_att(ncid,iVarId,'units',trim(unit_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

  if (present(vcal)) then ! add time calendar
    calendar_str = trim(vcal)
    ierr = nf90_put_att(ncid,iVarId,'calendar',trim(calendar_str))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
  end if

 END SUBROUTINE defvar


END MODULE write_simoutput
