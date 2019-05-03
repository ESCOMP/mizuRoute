module read_simoutput
USE nrtype
USE netcdf
implicit none
private
public::get_qDims
public::get_qMeta
public::getVarQsim
contains

 ! *********************************************************************
 ! new subroutine: read dimensions from runoff file
 ! *********************************************************************
 subroutine get_qDims(fname,           &  ! input: filename
                      vname_hruid,     &  ! input: name of coordinate dimension HRUid
                      vname_time,      &  ! input: name of coordinate dimension time
                      units_time,      &  ! output: time units
                      nTime,           &  ! output: number of time elements
                      nHRU,            &  ! output: number of HRUs in the runoff data file
                      ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_hruid  ! coordinate variable for HRUs
 character(*), intent(in)        :: vname_time   ! coordinate variable for time
 ! output variables
 character(*), intent(out)       :: units_time   ! time units
 integer(i4b), intent(out)       :: nTime        ! number of time elements
 integer(i4b), intent(out)       :: nHRU         ! number of HRUs
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 integer(i4b)                    :: idimID_time  ! dimension ID for time
 integer(i4b)                    :: idimID_hru   ! dimension ID for hru
 ! initialize error control
 ierr=0; message='get_qDims/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the ID of the time dimension
 ierr = nf90_inq_dimid(ncid, vname_time, idimID_time)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_time); return; endif

 ! get the length of the time dimension
 ierr = nf90_inquire_dimension(ncid, idimID_time, len=nTime)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the ID of the time variable
 ierr = nf90_inq_varid(ncid, trim(vname_time), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the time units
 ierr = nf90_get_att(ncid, ivarID, 'units', units_time)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the ID of the HRU dimension
 ierr = nf90_inq_dimid(ncid, vname_hruID, idimID_hru)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(vname_hruID); return; endif

 ! get the length of the HRU dimension
 ierr = nf90_inquire_dimension(ncid, idimID_hru, len=nHRU)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_qDims

 ! *********************************************************************
 ! new subroutine: read metadata from runoff file
 ! *********************************************************************
 subroutine get_qMeta(fname,           &  ! input: filename
                      vname_hruid,     &  ! input: name of coordinate dimension HRUid
                      qsimHRUid,       &  ! output: HRUid in the simulations
                      ierr, message)      ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_hruid  ! coordinate variable for HRUs
 ! output variables
 integer(i4b), intent(out)       :: qsimHRUid(:) ! HRUid in the simulations
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 ! initialize error control
 ierr=0; message='get_qMeta/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! get the variable ID for the HRUid
 !print*, 'trim(vname_hruid) = ', trim(vname_hruid)
 ierr = nf90_inq_varid(ncid, trim(vname_hruid), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the data for the HRUid
 !print*, 'size(qsimHRUid) = ', size(qsimHRUid)
 ierr = nf90_get_var(ncid, ivarID, qsimHRUid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//'; name='//trim(fname); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine get_qMeta


 ! *********************************************************************
 ! new subroutine: read runoff data
 ! *********************************************************************
 subroutine getVarQsim(fname,          &  ! input: filename
                       vname_time,     &  ! input: name of coordinate dimension time
                       vname_qsim,     &  ! input: name of runoff variable
                       iTime,          &  ! input: time index
                       dTime,          &  ! output: time
                       qsim_hru,       &  ! output: simulated runoff
                       ierr, message)     ! output: error control
 implicit none
 ! input variables
 character(*), intent(in)        :: fname        ! filename
 character(*), intent(in)        :: vname_time   ! name of coordinate dimension time
 character(*), intent(in)        :: vname_qsim   ! name of runoff variable
 integer(i4b), intent(in)        :: iTime        ! index of time element
 ! output variables
 real(dp), intent(out)           :: dTime        ! time
 real(dp), intent(out)           :: qsim_HRU(:)  ! simulated runoff for each HRU
 integer(i4b), intent(out)       :: ierr         ! error code
 character(*), intent(out)       :: message      ! error message
 ! local variables
 integer(i4b)                    :: ncid         ! NetCDF file ID
 integer(i4b)                    :: iVarID       ! variable ID
 real(dp),dimension(1)           :: tempTime     ! temporarary time vector (length=1)
 ! initialize error control
 ierr=0; message='getVarQsim/'

 ! open file for reading
 ierr = nf90_open(fname, nf90_nowrite, ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the variable ID for time
 ierr = nf90_inq_varid(ncid, trim(vname_time), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the time
 ierr = nf90_get_var(ncid, ivarID, tempTime, start=(/iTime/), count=(/1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
 dTime = tempTime(1)

 ! get the variable ID for the simulated runoff
 ierr = nf90_inq_varid(ncid, trim(vname_qsim), ivarID)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! get the simulated runoff data
 ierr = nf90_get_var(ncid, ivarID, qsim_HRU, start=(/1,iTime/), count=(/size(qsim_HRU),1/))
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 ! close the NetCDF file
 ierr = nf90_close(ncid)
 if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif

 end subroutine getVarQsim


end module read_simoutput
