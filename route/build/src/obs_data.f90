MODULE obs_data

! flow obs data class
! flow data should be stored in netcdf
!   dimensions:
!     time = nTime ;
!     site = nSite ;
!     stringLen = nChar ;
!   variables:
!     float streamflow(time, site)
!       units: m3/s
!     char  site(site, stringLen)
!     float time(time)
!       units: m3/s
!       calendar: regular
!
! usage:
! initialize gage obs data
!   gage_flow = gageObs(trim(ncname), varname_flow, varname_time, varname_site, dimname_time, dimname_site, ierr, cmessage)
!
! initialize reach water take data
!   qtake = waterTake(trim(ncname), varname_flow, varname_time, varname_reach, dimname_time, dimname_reach, ierr, cmessage)
!

  USE nrtype
  USE io_netcdf
  USE datetime_data, ONLY: datetime
  USE gageMeta_data, ONLY: gageMeta
  USE nr_utility_module, ONLY: match_index
  USE public_var, ONLY: integerMissing
  USE public_var, ONLY: clen=>strlen_gageSite

  implicit none

  private
  public:: gageObs
  public:: waterTake

  type, abstract :: obsData
    private
    integer(i4b)                   :: ncid = integerMissing ! netCDF ID
    character(strLen)              :: varName_obs
    character(strLen)              :: varName_time
    character(strLen)              :: varName_loc
    character(strLen)              :: dimName_time
    character(strLen)              :: dimName_loc
    integer(i4b)                   :: nLoc
    type(datetime),  allocatable   :: time(:)               ! datetime in the netcdf
    real(dp),        allocatable   :: obs(:,:)              ! selected observed data
    integer(i4b),    allocatable   :: link_index(:)         ! index of target location array (e.g., reachID)
    logical(lgt)                   :: fileOpen = .false.    ! flag to indicate history output netcdf is open

    CONTAINS

      procedure,  public :: openNC
      procedure,  public :: isFileOpen
      procedure,  public :: closeNC
      procedure,  public :: read_obs
      procedure,  public :: read_time
      procedure,  public :: time_ix => fn_get_time_ix
      procedure,  public :: link_Ix => fn_get_link_ix
      generic,    public :: get_datetime => fn_get_datetime, fn_get_datetime_vec, fn_get_datetime_scalar
      generic,    public :: get_obs => fn_get_obs, fn_get_obs_vec, fn_get_obs_tvec, fn_get_obs_svec, fn_get_obs_scalar
      procedure, private :: fn_get_link_ix
      procedure, private :: fn_get_time_ix
      procedure, private :: fn_get_datetime
      procedure, private :: fn_get_datetime_vec
      procedure, private :: fn_get_datetime_scalar
      procedure, private :: fn_get_obs
      procedure, private :: fn_get_obs_vec
      procedure, private :: fn_get_obs_tvec
      procedure, private :: fn_get_obs_svec
      procedure, private :: fn_get_obs_scalar

  end type obsData

  type, extends(obsData) :: gageObs
    character(25),   allocatable   :: site(:)  ! list of gauge site id
    CONTAINS
      procedure,  public :: read_site
      procedure,  public :: comp_link => sub_comp_link_site2reach
      generic,    public :: gage_id => fn_get_gageID, fn_get_gageID_vec, fn_get_gageID_scalar
      procedure,  public :: site_ix => fn_get_site_ix
      procedure, private :: fn_get_gageID
      procedure, private :: fn_get_gageID_vec
      procedure, private :: fn_get_gageID_scalar
  end type gageObs

  type, extends(obsData) :: waterTake
    integer(i4b),   allocatable   :: reach(:)  ! list of all the reach id
    CONTAINS
      procedure,  public :: read_reach
      procedure,  public :: comp_link => sub_comp_link_reach2reach
      generic,    public :: reach_id => fn_get_reachID, fn_get_reachID_vec, fn_get_reachID_scalar
      procedure,  public :: reach_ix => fn_get_reach_ix
      procedure, private :: fn_get_reachID
      procedure, private :: fn_get_reachID_vec
      procedure, private :: fn_get_reachID_scalar
  end type waterTake

  INTERFACE waterTake
    module procedure waterTake_constructor
  END INTERFACE waterTake

  INTERFACE gageObs
    module procedure gageObs_constructor
  END INTERFACE gageObs

  private:: fn_get_time_ix
  private:: fn_get_site_ix
  private:: fn_get_reach_ix
  private:: sub_comp_link_reach2reach
  private:: sub_comp_link_site2reach

  public :: get_time_ix
  public :: get_site_ix
  public :: get_reach_ix

  CONTAINS

    ! -----------------------------------------------------
    ! Instantiate gageObs
    ! -----------------------------------------------------
    FUNCTION gageObs_constructor(fileName, &
                                 obsVarName, &
                                 timeVarName, &
                                 locVarName, &
                                 timeDimName, &
                                 locDimName, &
                                 ierr, message) RESULT(instObs)
      implicit none
      type(gageObs)                        :: instObs
      character(*),           intent(in)   :: fileName
      character(*),           intent(in)   :: obsVarName
      character(*),           intent(in)   :: timeVarName
      character(*),           intent(in)   :: locVarName
      character(*),           intent(in)   :: timeDimName
      character(*),           intent(in)   :: locDimName
      integer(i4b),           intent(out)  :: ierr             ! error code
      character(*),           intent(out)  :: message          ! error message
      ! local variables
      character(strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='gageObs/'

      instObs%varName_obs  = obsVarName
      instObs%varName_time = timeVarName
      instObs%varName_loc  = locVarName
      instObs%dimName_time = timeDimName
      instObs%dimName_loc  = locDimName

      call instObs%openNC(fileName, ierr, cmessage) ! open gauge data nc and assing ncid and fileOpen = .true.
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call instObs%read_site(ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call instObs%read_time(ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      instObs%nLoc = size(instObs%site)

    END FUNCTION gageObs_constructor

    ! -----------------------------------------------------
    ! Instantiate gageObs
    ! -----------------------------------------------------
    FUNCTION waterTake_constructor(fileName, &
                                   obsVarName, &
                                   timeVarName, &
                                   locVarName, &
                                   timeDimName, &
                                   locDimName, &
                                   ierr, message) RESULT(instObs)
      implicit none
      type(waterTake)                      :: instObs
      character(*),           intent(in)   :: fileName
      character(*),           intent(in)   :: obsVarName
      character(*),           intent(in)   :: timeVarName
      character(*),           intent(in)   :: locVarName
      character(*),           intent(in)   :: timeDimName
      character(*),           intent(in)   :: locDimName
      integer(i4b),           intent(out)  :: ierr             ! error code
      character(*),           intent(out)  :: message          ! error message
      ! local variables
      character(strLen)                    :: cmessage         ! error message of downwind routine

      ierr=0; message='waterTake/'

      instObs%varName_obs  = obsVarName
      instObs%varName_time = timeVarName
      instObs%varName_loc  = locVarName
      instObs%dimName_time = timeDimName
      instObs%dimName_loc  = locDimName

      call instObs%openNC(fileName, ierr, cmessage) ! open gauge data nc and assing ncid and fileOpen = .true.
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call instObs%read_reach(ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call instObs%read_time(ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      instObs%nLoc = size(instObs%reach)

    END FUNCTION waterTake_constructor

    ! ---------------------------------
    ! open netCDF
    ! ---------------------------------
    SUBROUTINE openNC(this, fname, ierr, message)
      implicit none
      ! Argument variables
      class(obsData),   intent(inout)  :: this
      character(*),     intent(in)     :: fname
      integer(i4b),     intent(out)    :: ierr             ! error code
      character(*),     intent(out)    :: message          ! error message
      ! local variables
      character(strLen)                :: cmessage         ! error message of downwind routine

      ierr=0; message='openNC/'

      call open_nc(trim(fname), 'r', this%ncid, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      this%fileOpen = .true.

    END SUBROUTINE openNC

    ! ---------------------------------
    ! Check if file is open or not
    ! ---------------------------------
    logical(lgt) FUNCTION isFileOpen(this)
      implicit none
      class(obsData), intent(inout) :: this

      isFileOpen = this%fileOpen
    END FUNCTION

    ! ---------------------------------
    ! close netCDF
    ! ---------------------------------
    SUBROUTINE closeNC(this, ierr, message)
      implicit none
      class(obsData), intent(inout) :: this
      integer(i4b),     intent(out)    :: ierr             ! error code
      character(*),     intent(out)    :: message          ! error message
      ! local variables
      character(strLen)                :: cmessage         ! error message of downwind routine

      ierr=0; message='closeNC/'

      if (this%fileOpen) then
        call close_nc(this%ncid, ierr, cmessage)
        this%fileOpen = .false.
      endif
    END SUBROUTINE closeNC

    ! ---------------------------------
    ! get site ID
    ! ---------------------------------
    SUBROUTINE read_site(this, ierr, message, index_loc)

      implicit none
      ! Argument variables
      class(gageObs),          intent(inout) :: this
      integer(i4b),            intent(out)   :: ierr            ! error code
      character(*),            intent(out)   :: message         ! error message
      integer(i4b), optional,  intent(in)    :: index_loc(:)
      ! local variables
      integer(i4b)                           :: nSite           ! total number of sites
      integer(i4b)                           :: nSite_read      ! number of sites read in
      integer(i4b)                           :: ix
      character(clen), allocatable           :: array_temp(:)
      character(strLen)                      :: cmessage        ! error message of downwind routine

      ierr=0; message='read_site/'

      call get_nc_dim_len(this%ncid, this%dimName_loc, nSite, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (present(index_loc)) then
        nSite_read = size(index_loc)

        allocate(array_temp(nSite), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

        allocate(this%site(nSite_read), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%site]'; return; endif

        call get_nc(this%ncid, this%varName_loc, array_temp, [1,1], [clen, nSite], ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        do ix=1,nSite_read
          this%site(ix) = array_temp(index_loc(ix))
        end do
      else
        allocate(this%site(nSite), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%site]'; return; endif

        call get_nc(this%ncid, this%varName_loc, this%site, [1,1], [clen, nSite], ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

    END SUBROUTINE read_site

    ! ---------------------------------
    ! get reach ID
    ! ---------------------------------
    SUBROUTINE read_reach(this, ierr, message, index_loc)

      implicit none
      ! Argument variables
      class(waterTake),        intent(inout) :: this
      integer(i4b),            intent(out)   :: ierr            ! error code
      character(*),            intent(out)   :: message         ! error message
      integer(i4b), optional,  intent(in)    :: index_loc(:)
      ! local variables
      integer(i4b)                           :: nSite           ! total number of sites
      integer(i4b)                           :: nSite_read      ! number of sites read in
      integer(i4b)                           :: ix
      integer(i4b), allocatable              :: array_temp(:)
      character(strLen)                      :: cmessage        ! error message of downwind routine

      ierr=0; message='read_reach/'

      call get_nc_dim_len(this%ncid, this%dimName_loc, nSite, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (present(index_loc)) then
        nSite_read = size(index_loc)

        allocate(array_temp(nSite), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

        allocate(this%reach(nSite_read), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%reach]'; return; endif

        call get_nc(this%ncid, this%varName_loc, array_temp, 1, nSite, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        do ix=1,nSite_read
          this%reach(ix) = array_temp(index_loc(ix))
        end do
      else
        allocate(this%reach(nSite), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%reach]'; return; endif

        call get_nc(this%ncid, this%varName_loc, this%reach, 1, nSite, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

    END SUBROUTINE read_reach

    ! ---------------------------------
    ! get time
    ! ---------------------------------
    SUBROUTINE read_time(this, ierr, message, index_time)

      implicit none
      ! Argument variables
      class(obsData),         intent(inout) :: this
      integer(i4b),           intent(out)   :: ierr            ! error code
      character(*),           intent(out)   :: message         ! error message
      integer(i4b), optional, intent(in)    :: index_time(:)
      ! local variables
      type(datetime)                        :: datetime_start
      real(dp), allocatable                 :: timeVar_temp(:)
      real(dp)                              :: timeVar
      real(dp)                              :: convTime2Secs
      integer(i4b)                          :: nTime           ! total number of sites
      integer(i4b)                          :: nTime_read      ! number of sites read in
      integer(i4b)                          :: ix
      character(100)                        :: timestr_start
      character(30)                         :: gage_calendar
      character(30)                         :: t_unit
      character(len=strLen)                 :: cmessage        ! error message of downwind routine

      ierr=0; message='read_time/'

      call get_nc_dim_len(this%ncid, this%dimName_time, nTime, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call get_var_attr(this%ncid, this%varName_time, 'units', timestr_start, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call get_var_attr(this%ncid, this%varName_time, 'calendar', gage_calendar, ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      call datetime_start%str2datetime(timestr_start, ierr, cmessage)

      t_unit = trim(timestr_start(1:index(timestr_start,' ')))
      select case( trim(t_unit) )
        case('seconds','second','sec','s'); convTime2Secs=1._dp
        case('minutes','minute','min','m'); convTime2Secs=24._dp
        case('hours'  ,'hour'  ,'hr' ,'h'); convTime2Secs=1440._dp
        case('days'   ,'day'   ,'d');       convTime2Secs=86400._dp
        case default
          ierr=20; message=trim(message)//'<t_unit>= '//trim(t_unit)//': <t_unit> must be seconds, minutes, hours or days.'; return
       end select

      if (present(index_time)) then
        nTime_read = size(index_time)

        allocate(timeVar_temp(nTime), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [timeVar_temp]'; return; endif

        allocate(this%time(nTime_read), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%time]'; return; endif

        call get_nc(this%ncid, this%varName_time, timeVar_temp, 1, nTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        do ix=1,nTime_read
          timeVar = timeVar_temp(index_time(ix))*convTime2Secs
          this%time(ix) = datetime_start%add_sec(timeVar, gage_calendar, ierr, cmessage)
        end do
      else
        allocate(timeVar_temp(nTime), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%time]'; return; endif

        allocate(this%time(nTime), stat=ierr, errmsg=cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%time]'; return; endif

        call get_nc(this%ncid, this%varName_time, timeVar_temp, 1, nTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! construct datetime data
        do ix=1,nTime
          timeVar = timeVar_temp(ix)*convTime2Secs
          this%time(ix) = datetime_start%add_sec(timeVar, gage_calendar, ierr, cmessage)
        end do

      end if

    END SUBROUTINE read_time

    ! ---------------------------------
    ! get obs data
    ! ---------------------------------
    SUBROUTINE read_obs(this, ierr, message, index_loc, index_time)

      implicit none
      ! Argument variables
      class(obsData),          intent(inout) :: this
      integer(i4b),            intent(out)   :: ierr             ! error code
      character(*),            intent(out)   :: message          ! error message
      integer(i4b), optional,  intent(in)    :: index_loc
      integer(i4b), optional,  intent(in)    :: index_time
      ! local variables
      real(dp),    allocatable               :: array_temp(:,:)
      integer(i4b)                           :: nLoc_read        ! number of sites read
      integer(i4b)                           :: nTime_read       ! number of time read
      integer(i4b)                           :: nTime            ! number of time
      integer(i4b)                           :: nLoc             ! number of site
      integer(i4b)                           :: ix_loc, ix_time  ! starting reading index
      integer(i4b)                           :: ix,jx
      character(len=strLen)                  :: cmessage         ! error message of downwind routine

      ierr=0; message='read_obs/'

      if (present(index_loc)) then
        nLoc_read = 1
        ix_loc = index_loc
      else
        call get_nc_dim_len(this%ncid, this%dimName_loc, nLoc, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        nLoc_read = nLoc
        ix_loc = 1
      end if
      if (present(index_time)) then
        nTime_read = 1
        ix_time = index_time
      else
        call get_nc_dim_len(this%ncid, this%dimName_time, nTime, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
        nTime_read = nTime
        ix_time = 1
      end if

      ! read full obs data
      allocate(array_temp(nLoc_read, nTime_read), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [array_temp]'; return; endif

      call get_nc(this%ncid, this%varName_obs, array_temp, [ix_loc, ix_time], [nLoc_read, nTime_read], ierr, cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      if (allocated(this%obs)) then
        deallocate(this%obs)
      end if
      allocate(this%obs(nTime_read, nLoc_read), stat=ierr, errmsg=cmessage)
      if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [this%obs]'; return; endif
      do ix=1,nLoc_read
        do jx=1,nTime_read
          this%obs(jx, ix) = array_temp(ix, jx)
        end do
      end do

    END SUBROUTINE read_obs

    ! -----------------------------------------------------
    ! Wrapper FUNCTIONs: get array index
    ! -----------------------------------------------------
    FUNCTION fn_get_time_ix(this, datetime_in) result(ix)
      implicit none
      class(obsData), intent(in) :: this
      type(datetime), intent(in) :: datetime_in
      integer(i4b)               :: ix
      ix = get_time_ix(this%time, datetime_in)
    END FUNCTION fn_get_time_ix

    FUNCTION fn_get_site_ix(this, site_in) result(ix)
      implicit none
      class(gageObs), intent(in) :: this
      character(*), intent(in)   :: site_in
      integer(i4b)               :: ix
      ix = get_site_ix(this%site, site_in)
    END FUNCTION fn_get_site_ix

    FUNCTION fn_get_reach_ix(this, reach_in) result(ix)
      implicit none
      class(waterTake), intent(in) :: this
      integer(i4b), intent(in)    :: reach_in
      integer(i4b)                :: ix
      ix = get_reach_ix(this%reach, reach_in)
    END FUNCTION fn_get_reach_ix

    ! -----------------------------------------------------
    ! FUNCTIONs: get gage ID
    ! -----------------------------------------------------
    pure FUNCTION fn_get_link_ix(this) result(ix)
      implicit none
      class(obsData), intent(in) :: this
      integer(i4b), allocatable  :: ix(:)
      allocate(ix(size(this%link_index)))
      ix = this%link_index
    END FUNCTION fn_get_link_ix

    ! -----------------------------------------------------
    ! FUNCTIONs: get gage ID
    ! -----------------------------------------------------
    pure FUNCTION fn_get_gageID(this) result(gageID)
      implicit none
      class(gageObs), intent(in) :: this
      character(clen), allocatable  :: gageID(:)
      allocate(gageID(size(this%site)))
      gageID = this%site
    END FUNCTION fn_get_gageID

    pure FUNCTION fn_get_gageID_vec(this, ix) result(gageID)
      implicit none
      class(gageObs), intent(in) :: this
      integer(i4b), intent(in)   :: ix(:)
      character(clen), allocatable  :: gageID(:)
      allocate(gageID(size(ix)))
      gageID = this%site(ix)
    END FUNCTION fn_get_gageID_vec

    pure FUNCTION fn_get_gageID_scalar(this, ix) result(gageID)
      implicit none
      class(gageObs), intent(in) :: this
      integer(i4b), intent(in)   :: ix
      character(clen)            :: gageID
      gageID = this%site(ix)
    END FUNCTION fn_get_gageID_scalar

    ! -----------------------------------------------------
    ! FUNCTIONs: get reach ID
    ! -----------------------------------------------------
    pure FUNCTION fn_get_reachID(this) result(reachID)
      implicit none
      class(waterTake), intent(in) :: this
      integer(i4b), allocatable    :: reachID(:)
      allocate(reachID(size(this%reach)))
      reachID = this%reach
    END FUNCTION fn_get_reachID

    pure FUNCTION fn_get_reachID_vec(this, ix) result(reachID)
      implicit none
      class(waterTake), intent(in) :: this
      integer(i4b), intent(in)     :: ix(:)
      integer(i4b), allocatable    :: reachID(:)
      allocate(reachID(size(ix)))
      reachID = this%reach(ix)
    END FUNCTION fn_get_reachID_vec

    pure FUNCTION fn_get_reachID_scalar(this, ix) result(reachID)
      implicit none
      class(waterTake), intent(in) :: this
      integer(i4b), intent(in)     :: ix
      integer(i4b)                 :: reachID
      reachID = this%reach(ix)
    END FUNCTION fn_get_reachID_scalar

    ! -----------------------------------------------------
    ! FUNCTIONs: get time
    ! -----------------------------------------------------
    FUNCTION fn_get_datetime(this) result(out_datetime)
      implicit none
      class(obsData),  intent(in) :: this
      type(datetime), allocatable :: out_datetime(:)
      integer(i4b)                :: ix
      allocate(out_datetime(size(this%time)))
      do ix=1,size(this%time)
        out_datetime(ix) = this%time(ix)
      end do
    END FUNCTION fn_get_datetime

    FUNCTION fn_get_datetime_vec(this, ix) result(out_datetime)
      implicit none
      class(obsData), intent(in)  :: this
      integer(i4b), intent(in)    :: ix(:)
      type(datetime), allocatable :: out_datetime(:)
      integer(i4b)                :: jx
      allocate(out_datetime(size(ix)))
      do jx=1,size(ix)
        out_datetime(jx) = this%time(ix(jx))
      end do
    END FUNCTION fn_get_datetime_vec

    FUNCTION fn_get_datetime_scalar(this, ix) result(out_datetime)
      implicit none
      class(obsData), intent(in) :: this
      integer(i4b),   intent(in) :: ix
      type(datetime)           :: out_datetime
      out_datetime = this%time(ix)
    END FUNCTION fn_get_datetime_scalar

    ! -----------------------------------------------------
    ! FUNCTIONs: get flow data
    ! -----------------------------------------------------
    pure FUNCTION fn_get_obs(this) result(flow)
      implicit none
      class(obsData), intent(in) :: this
      real(dp), allocatable      :: flow(:,:)
      allocate(flow(size(this%obs,1), size(this%obs,2)))
      flow = this%obs
    END FUNCTION fn_get_obs

    pure FUNCTION fn_get_obs_vec(this, tix, six) result(flow)
      implicit none
      class(obsData), intent(in) :: this
      integer(i4b),   intent(in) :: tix(:)
      integer(i4b),   intent(in) :: six(:)
      real(dp), allocatable     :: flow(:,:)
      allocate(flow(size(tix),size(six)))
      flow = this%obs(tix, six)
    END FUNCTION fn_get_obs_vec

    pure FUNCTION fn_get_obs_tvec(this, tix, six) result(flow)
      implicit none
      class(obsData),intent(in) :: this
      integer(i4b),  intent(in) :: tix(:)
      integer(i4b),  intent(in) :: six
      real(dp), allocatable     :: flow(:)
      allocate(flow(size(tix)))
      flow = this%obs(tix, six)
    END FUNCTION fn_get_obs_tvec

    pure FUNCTION fn_get_obs_svec(this, tix, six) result(flow)
      implicit none
      class(obsData), intent(in) :: this
      integer(i4b),   intent(in) :: tix
      integer(i4b),   intent(in) :: six(:)
      real(dp), allocatable      :: flow(:)
      allocate(flow(size(six)))
      flow = this%obs(tix, six)
    END FUNCTION fn_get_obs_svec

    pure FUNCTION fn_get_obs_scalar(this, tix, six) result(flow)
      implicit none
      class(obsData), intent(in) :: this
      integer(i4b),  intent(in) :: tix
      integer(i4b),  intent(in) :: six
      real(dp)                  :: flow
      flow = this%obs(tix,six)
    END FUNCTION fn_get_obs_scalar

    ! -----------------------------------------------------
    ! subroutine: index of reach ID array given reach iD
    ! -----------------------------------------------------
    SUBROUTINE sub_comp_link_reach2reach(this, reachID_in)
      implicit none
      ! Argument variables
      class(waterTake),       intent(inout)  :: this
      integer(i4b),           intent(in)     :: reachID_in(:)

      allocate(this%link_index(this%nLoc))

      ! Find index of matching reachID in reachID_in array
      this%link_index = match_index(reachID_in, this%reach, missingValue=integerMissing)
    END SUBROUTINE sub_comp_link_reach2reach

    ! -----------------------------------------------------
    ! subroutine: index of reach ID array given site iD
    ! -----------------------------------------------------
    SUBROUTINE sub_comp_link_site2reach(this, reachID_in, gageMeta_in)
      implicit none
      ! Argument variables
      class(gageObs),         intent(inout)  :: this
      integer(i4b),           intent(in)     :: reachID_in(:)
      type(gageMeta),         intent(in)     :: gageMeta_in
      ! local variables
      integer(i4b), allocatable              :: reachID_gage(:)
      integer(i4b)                           :: ix, jx

      allocate(this%link_index(this%nLoc), reachID_gage(this%nLoc))

      ! assuming nLoc and nGage are small
      do ix = 1, this%nLoc
        do jx = 1, gageMeta_in%gage_num()
          if (trim(this%site(ix)) == trim(gageMeta_in%gage_id(jx))) then
            reachID_gage(ix) = gageMeta_in%reach_id(jx)
            exit
          end if
        end do
      end do

      ! Find index of matching reachID in reachID_in array
      this%link_index = match_index(reachID_in, reachID_gage, missingValue=integerMissing)
    END SUBROUTINE sub_comp_link_site2reach

    ! -----------------------------------------------------
    ! private functions/subroutines
    ! -----------------------------------------------------
    FUNCTION get_time_ix(datetime_array, datetime_in) result(tix)
      implicit none
      type(datetime), intent(in) :: datetime_array(:)
      type(datetime), intent(in) :: datetime_in
      integer(i4b)   :: tix
      integer(i4b)   :: ix
      tix=integerMissing
      do ix = 1, size(datetime_array)
        if (datetime_array(ix) == datetime_in) then
          tix = ix
          exit
        end if
      end do
    END FUNCTION get_time_ix

    FUNCTION get_site_ix(site_array, site_in) result(six)
      implicit none
      character(*), intent(in) :: site_array(:)
      character(*), intent(in) :: site_in
      integer(i4b)   :: six
      integer(i4b)   :: ix
      six=integerMissing
      do ix = 1, size(site_array)
        if (trim(site_array(ix)) == trim(site_in)) then
          six = ix
          exit
        end if
      end do
    END FUNCTION get_site_ix

    FUNCTION get_reach_ix(reach_array, reach_in) result(rix)
      implicit none
      integer(i4b), intent(in) :: reach_array(:)
      integer(i4b), intent(in) :: reach_in
      integer(i4b) :: rix
      integer(i4b) :: ix
      rix=integerMissing
      do ix = 1, size(reach_array)
        if (reach_array(ix) == reach_in) then
          rix = ix
          exit
        end if
      end do
    END FUNCTION get_reach_ix

END MODULE obs_data
