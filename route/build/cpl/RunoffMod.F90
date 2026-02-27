MODULE RunoffMod

! DESCRIPTION:
! Module containing data structures for coupler runoff data

  USE nrtype,     ONLY : r8 => dp, cs => strLen
  USE public_var, ONLY : iulog
  USE public_var, ONLY : integerMissing
  USE shr_sys_mod,ONLY : shr_sys_abort, shr_sys_flush

  implicit none

  private

  type, public :: rof_control
    ! tracers
    integer           :: ntracers = integerMissing   ! number of tracers
    character(len=cs), allocatable :: tracer_names(:) ! tracer names
    integer           :: nt_liq                      ! index of liquid tracer in tracer_names
    integer           :: nt_ice                      ! index of ice tracer in tracer_names
    logical           :: rof_from_glc                ! if true, will receive liq and ice runoff from glc
    ! domain and decomposition info
    integer           :: begr,endr        ! local start/stop indices
    integer           :: lnumr            ! local number of catchments
    integer           :: numr             ! gdc global number of catchments
    integer , pointer :: gindex(:)        ! global index consistent with map file
    real(r8), pointer :: area(:)          ! area of hru (AKA catchment) [m2]
    integer , pointer :: mask(:)          ! reach category. 1=non outle, 2=interior basin outlet, 3=ocean outlet
    ! import variables
    real(r8), pointer :: qsur(:,:)        ! coupler importing surface forcing [mm/s]
    real(r8), pointer :: qsub(:,:)        ! coupler importing subsurface forcing [mm/s]
    real(r8), pointer :: qgwl(:,:)        ! coupler importing glacier/wetland/lake forcing [mm/s]
    real(r8), pointer :: qirrig(:)        ! coupler importing irrigation [mm/s] - negative
    real(r8), pointer :: qirrig_actual(:) ! actual water take from river reach [mm/s]
    ! export variables
    real(r8), pointer :: direct(:,:)      ! coupler return direct flow to ocean [mm/s]
    real(r8), pointer :: discharge(:,:)   ! coupler exporting river discharge [m3/s]
    real(r8), pointer :: volr(:)          ! coupler exporting river storage per unit HRU area (m)
    real(r8), pointer :: flood(:)         ! coupler exporting flood water sent back to clm [m3/s]
  CONTAINS
    procedure, public  :: init
    procedure, public  :: init_tracer_names
  end type rof_control

CONTAINS

  SUBROUTINE init(this, begr, endr, numr)
    implicit none
    class(rof_control)  :: this
    integer, intent(in) :: begr, endr, numr
    integer :: nt
    integer :: ierr

    nt = this%ntracers
    allocate(this%gindex(begr:endr),          &
             this%area(begr:endr),            &
             this%mask(begr:endr),            &
             this%qsur(begr:endr,nt),         &
             this%qsub(begr:endr,nt),         &
             this%qgwl(begr:endr,nt),         &
             this%qirrig(begr:endr),          &
             this%qirrig_actual(begr:endr),   &
             this%discharge(begr:endr,nt),    &
             this%direct(begr:endr,nt),       &
             this%volr(begr:endr),            &
             this%flood(begr:endr),           &
             stat=ierr)
    if (ierr/=0) then
      write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
      call shr_sys_flush(iulog)
      call shr_sys_abort
    end if

    this%gindex(:)       = integerMissing
    this%mask(:)         = integerMissing
    this%area(:)         = 0._r8
    this%qirrig(:)       = 0._r8
    this%qirrig_actual(:)= 0._r8
    this%qsur(:,:)       = 0._r8
    this%qsub(:,:)       = 0._r8
    this%qgwl(:,:)       = 0._r8
    this%discharge(:,:)  = 0._r8
    this%direct(:,:)     = 0._r8
    this%volr(:)         = 0._r8
    this%flood(:)        = 0._r8

  END SUBROUTINE init

  SUBROUTINE init_tracer_names(this, tracer_names)
    implicit none
    ! Argument variables
    class(rof_control)            :: this
    character(len=cs), intent(in) :: tracer_names(:)    ! string array of tracer names
    ! Local variables
    integer :: nt
    character(len=*),parameter :: subname = '(RunoffMod: init_tracer_names)'

    ! Determine number of tracers and array of tracer names
    this%ntracers = len(tracer_names)
    allocate(this%tracer_names(this%ntracers))
    do nt = 1,this%ntracers
      this%tracer_names(nt) = tracer_names(nt)
    end do
    ! Set tracers
    this%nt_liq = 0
    this%nt_ice = 0
    do nt = 1,this%ntracers
       if (trim(this%tracer_names(nt)) == 'LIQ') this%nt_liq = nt
       if (trim(this%tracer_names(nt)) == 'ICE') this%nt_ice = nt
    enddo
    if (this%nt_liq == 0 .or. this%nt_ice == 0) then
       write(iulog,*) trim(subname),': ERROR in tracers LIQ ICE ',this%nt_liq,this%nt_ice,this%tracer_names(:)
       call shr_sys_abort()
    endif
  END SUBROUTINE init_tracer_names

END MODULE RunoffMod
