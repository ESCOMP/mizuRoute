module get_runoff

USE nrtype
USE dataTypes, ONLY: infileinfo

USE io_netcdf, only:get_nc
USE io_netcdf, only:get_var_attr_real
USE io_netcdf, only:get_var_attr_char

USE io_netcdf, only:get_nc_dim_len

implicit none

private

public::get_hru_runoff

contains

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine get_hru_runoff(ierr, message)     ! output: error control

 ! populate runoff_data with runoff values at LSM domain and at iTime step

  ! subroutines
  USE process_time_module, ONLY: process_time  ! process time information
  USE read_runoff, only:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   only:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   only:sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)
  ! external subroutines
  USE ascii_util_module, ONLY:file_open        ! open file (performs a few checks as well)
  USE ascii_util_module, ONLY:get_vlines       ! get a list of character strings from non-comment lines


  ! shared data
  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE public_var,  only:vname_time              ! public var time varibale name
  USE public_var,  only:time_units              ! time units (seconds, hours, or days) from control file
  USE public_var,  only: simStart               ! date string defining the start of the simulation
  USE public_var,  only: simEnd                 ! date string defining the end of the simulation
  USE globalData,  only:iTime
  USE globalData,  only:nHRU
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data
  USE globalData,  only: timeVar                ! time variables (unit given by runoff data)
  USE globalData,  only: convTime2Days          ! conversion multipliers for time unit of runoff input to day
  USE public_var,  only: calendar               ! calendar name
  USE globalData,  only: refJulday               ! julian day: reference
  USE globalData,  only: startJulday             ! julian day: start of routing simulation
  USE globalData,  only: endJulday               ! julian day: end of routing simulation

  implicit none
  ! input variables: none
  ! output variables
  integer(i4b), intent(out)     :: ierr                                 ! error code
  character(*), intent(out)     :: message                              ! error message
  ! local variables
  real(dp)    , allocatable     :: basinRunoff(:)                       ! basin runoff (m/s)
  real(dp)    , allocatable     :: basinEvapo(:)                        ! basin runoff (m/s)
  real(dp)    , allocatable     :: basinPrecip(:)                       ! basin runoff (m/s)
  character(len=strLen)         :: runoff = 'runoff'                    ! flag in case the flux is runoff
  character(len=strLen)         :: evaporation = 'evaporation'          ! flag in case the flux is evaporation
  character(len=strLen)         :: precipitation = 'precipitation'      ! flag in case the flux is precipitation
  character(len=strLen)         :: cmessage                             ! error message from subroutine
  integer(i4b)                  :: size_fname_qsim                      ! error code
  character(len=strLen)         :: infile                               ! input filename
  integer(i4b)                  :: unt                                  ! file unit (free unit output from file_open)
  character(len=strLen),allocatable :: dataLines(:)                     ! vector of lines of information (non-comment lines)
  integer(i4b)                      :: iFile                            ! counter for forcing files
  integer(i4b)                      :: nFile                            ! number of forcing files in forcing file list
  character(len=strLen)             :: filenameData                     ! name of forcing datafile
  type(infileinfo), allocatable     ::  infileinfo(:)                   ! the file
  integer(i4b), parameter           ::  nTime = 4662                    ! hard coded for now

  ! initialize error control
  ierr=0; message='get_hru_runoff/'

  fname_qsim = 'NetCDF_input.txt' ! just chenge the file name locally here to check the model simulation
  print*, 'print default fname with input in get_basin_runoff',fname_qsim

  print*, 'iTime is in get_basin_runoff:',iTime
  ! checking if fname_qsim is .txt or .nc
  size_fname_qsim = len(trim(fname_qsim))
  print*, size_fname_qsim
  print*, 'the extension of the file',fname_qsim(size_fname_qsim-3:size_fname_qsim)
  print*, 'the fname_qsim', fname_qsim
  if (fname_qsim (size_fname_qsim-3:size_fname_qsim).eq.'.txt')then
   print*, '*.txt'

    ! ------------------------------------------------------------------------------------------------------------------
    ! (1) read from the list of forcing files and populate the start and end time steps in Julian days
    ! ------------------------------------------------------------------------------------------------------------------
    ! build filename for forcing file list
    infile = trim(input_dir)//trim(fname_qsim)

    ! open file
    call file_open(trim(infile),unt,ierr,cmessage)  ! question... what is the unt?
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

    ! get a list of character strings from non-comment lines
    call get_vlines(unt,dataLines,ierr,cmessage)
    if(ierr/=0)then; ierr=20; message=trim(message)//trim(cmessage); return; end if
    nFile = size(dataLines)
    print*, 'number of lines in the text file', nFile

    ! allocate space for forcing information
    if(allocated(infileinfo)) deallocate(infileinfo)
    allocate(infileinfo(nFile), stat=ierr)
    if(ierr/=0)then; ierr=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; end if

    ! poputate the forcingInfo structure with filenames
    do iFile=1,nFile
     ! split the line into "words" (expect one word: the file describing forcing data for that index)
     read(dataLines(iFile),*,iostat=ierr) filenameData
     if(ierr/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infile)//']'; return; end if
     ! set forcing file name attribute
     infileinfo(iFile)%infilename = trim(filenameData)
     !
     call get_nc(trim(input_dir)//trim(fname_qsim), vname_time, timeVar, 1, nTime, ierr, cmessage)

     ! get the time multiplier needed to convert time to units of days
     select case( trim( time_units(1:index(time_units,' ')) ) )
      case('seconds'); convTime2Days=86400._dp
      case('hours');   convTime2Days=24._dp
      case('days');    convTime2Days=1._dp
      case default;    ierr=20; message=trim(message)//'unable to identify time units'; return
     end select
     ! vname_time = vname_time / convTime2Days

     print*, time_units
     ! extract time information from the control information
     call process_time(time_units,    calendar, refJulday,   ierr, cmessage) ! the time unit is coming from the NetCDF for one NetCDF file
     if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [refJulday]'; return; endif
     ! similarly the endref should be done here....
     ! refJuldayend = refJulday + vname_time (end)

     !! startJulyday and endJulyday can be move to earlier part as they are used only once...
     call process_time(trim(simStart),calendar, startJulday, ierr, cmessage) ! start of the simulation from the
     if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [startJulday]'; return; endif
     call process_time(trim(simEnd),  calendar, endJulday,   ierr, cmessage)
     if(ierr/=0) then; message=trim(message)//trim(cmessage)//' [endJulday]'; return; endif

     print*, refJulday, startJulday, endJulday

  ! check that the dates are aligned
  if(endJulday<startJulday) then; ierr=20; message=trim(message)//'simulation end is before simulation start'; return; endif

    end do  ! (looping through files)
    print*, infileinfo(1)%infilename
    print*, infileinfo(2)%infilename

    ! close ascii file
    close(unit=unt,iostat=ierr); if(ierr/=0)then;message=trim(message)//'problem closing forcing file list'; return; end if

  endif
  if (fname_qsim (size_fname_qsim-2:size_fname_qsim).eq.'.nc')then
   print*, '*.nc'
  endif
  stop

  ! get the simulated runoff for the current time step - runoff_data%qsim(:), %qsim2D(:,:), easim(:), easim2d(:,:), precip(:) and precip2d(:)
  call read_runoff_data(trim(input_dir)//trim(fname_qsim), & ! input: filename
                        iTime,                             & ! input: time index
                        runoff_data,                       & ! inout: runoff data structure
                        ierr, cmessage)                      ! output: error control

  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinRunoff (local array)
  allocate(basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinRunoff
  if ( allocated(runoff_data%basinRunoff) ) then
    deallocate(runoff_data%basinRunoff, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinRunoff(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinEvapo (local array)
  allocate(basinEvapo(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinEvapo
  if ( allocated(runoff_data%basinEvapo) ) then
    deallocate(runoff_data%basinEvapo, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinEvapo(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate basinPrecip (local array)
  allocate(basinPrecip(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! initialize runoff_data%basinPrecip
  if ( allocated(runoff_data%basinPrecip) ) then
    deallocate(runoff_data%basinPrecip, stat=ierr)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  end if
  allocate(runoff_data%basinPrecip(nHRU), stat=ierr)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif


  ! Get river network HRU runoff into runoff_data data structure
  if (is_remap) then ! remap LSM simulated runoff to the HRUs in the river network

   call remap_runoff(runoff_data, remap_data, runoff,        basinRunoff, ierr, cmessage) ! if flage runoff then it uses qsim or qsim2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call remap_runoff(runoff_data, remap_data, evaporation,   basinEvapo, ierr, cmessage) ! if flag evaporation then it uses easim or easim2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
   call remap_runoff(runoff_data, remap_data, Precipitation, basinPrecip, ierr, cmessage) ! if flag precipitation then it uses precip or precip2d
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  else ! runoff is already remapped to river network HRUs

   call sort_runoff(runoff_data, runoff, basinRunoff, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call sort_runoff(runoff_data, evaporation, basinEvapo,  ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call sort_runoff(runoff_data, precipitation, basinPrecip, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  end if

  runoff_data%basinRunoff = basinRunoff
  runoff_data%basinEvapo = basinEvapo
  runoff_data%basinPrecip = basinPrecip

 end subroutine get_hru_runoff

end module get_runoff
