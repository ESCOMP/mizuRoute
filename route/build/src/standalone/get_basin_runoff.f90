module get_runoff

USE nrtype
USE dataTypes, ONLY: infileinfo

implicit none

private

public::get_hru_runoff

contains

 ! *********************************************************************
 ! public subroutine: read runoff data
 ! *********************************************************************
 subroutine get_hru_runoff(ierr, message)     ! output: error control

 ! populate runoff_data with runoff values at LSM domain and at iTime step

  ! shared data
  USE public_var,  only:input_dir               ! directory containing input data
  USE public_var,  only:fname_qsim              ! simulated runoff netCDF name
  USE public_var,  only:is_remap                ! logical whether or not runnoff needs to be mapped to river network HRU
  USE globalData,  only:iTime
  USE globalData,  only:nHRU
  USE globalData,  only:runoff_data             ! data structure to hru runoff data
  USE globalData,  only:remap_data              ! data structure to remap data
  ! subroutines
  USE read_runoff, only:read_runoff_data        ! read runoff value into runoff_data data strucuture
  USE remapping,   only:remap_runoff            ! mapping HM runoff to river network HRU runoff (HM_HRU /= RN_HRU)
  USE remapping,   only:sort_runoff             ! mapping HM runoff to river network HRU runoff (HM_HRU == RN_HRU)
  ! external subroutines
  USE ascii_util_module, ONLY:file_open        ! open file (performs a few checks as well)
  USE ascii_util_module, ONLY:get_vlines       ! get a list of character strings from non-comment lines

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
  type(infileinfo), allocatable ::  infileinfo(:)                       ! the file

  ! initialize error control
  ierr=0; message='get_hru_runoff/'


  print*, 'iTime is:',iTime
  ! checking if fname_qsim is .txt or .nc
  size_fname_qsim = len(trim(fname_qsim))
  print*, size_fname_qsim
  print*, fname_qsim(size_fname_qsim-3:size_fname_qsim)
  print*, fname_qsim
  if (fname_qsim (size_fname_qsim-3:size_fname_qsim).eq.'.txt')then
   print*, '*.txt'

    ! ------------------------------------------------------------------------------------------------------------------
    ! (1) read from the list of forcing files
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
    end do  ! (looping through files)

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
