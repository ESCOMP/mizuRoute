MODULE RunoffMod

! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing data structures for coupler runoff data

  USE nrtype,     ONLY : r8 => dp
  USE public_var, ONLY : iulog, MaxPosVal
  USE rtmVar,     ONLY : nt_rtm
  USE shr_sys_mod,ONLY : shr_sys_abort

! !PUBLIC TYPES:
  implicit none

  private

  type, public :: runoff_flow
     integer           :: begr,endr        ! local start/stop indices
     integer           :: lnumr            ! local number of catchments
     integer           :: numr             ! gdc global number of catchments

     integer , pointer :: gindex(:)        ! global index consistent with map file

     real(r8), pointer :: area(:)          ! area of catchment
     integer , pointer :: mask(:)          ! reach category. 1=non outle, 2=interior basin outlet, 3=ocean outlet

     real(r8), pointer :: qsur(:,:)        ! coupler importing surface forcing [m3/s]
     real(r8), pointer :: qsub(:,:)        ! coupler importing subsurface forcing [m3/s]
     real(r8), pointer :: qgwl(:,:)        ! coupler importing glacier/wetland/lake forcing [m3/s]
     real(r8), pointer :: qirrig(:)        ! coupler importing irrigation [m3/s]
     real(r8), pointer :: qirrig_actual(:) ! minimum of irrigation and available main channel storage

     real(r8), pointer :: discharge(:,:)   ! coupler exporting river discharge [m3/s]
     real(r8), pointer :: volr(:)          ! coupler exporting river storage (m3)
     real(r8), pointer :: flood(:)         ! coupler exporting flood water sent back to clm [m3/s]
  end type runoff_flow

  type(runoff_flow), public :: rtmCTL

  public :: RunoffInit

CONTAINS

  SUBROUTINE RunoffInit(begr, endr, numr)

    integer, intent(in) :: begr, endr, numr

    integer :: ierr

    allocate(rtmCTL%gindex(begr:endr),            &
             rtmCTL%area(begr:endr),              &
             rtmCTL%mask(begr:endr),              &
             rtmCTL%qsur(begr:endr,nt_rtm),       &
             rtmCTL%qsub(begr:endr,nt_rtm),       &
             rtmCTL%qgwl(begr:endr,nt_rtm),       &
             rtmCTL%qirrig(begr:endr),            &
             rtmCTL%qirrig_actual(begr:endr),     &
             rtmCTL%discharge(begr:endr,nt_rtm),  &
             rtmCTL%volr(begr:endr),              &
             rtmCTL%flood(begr:endr),             &
             stat=ierr)
    if (ierr /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff local arrays'
       call shr_sys_abort
    end if

    rtmCTL%qirrig(:)       = 0._r8
    rtmCTL%qirrig_actual(:)= 0._r8
    rtmCTL%qsur(:,:)       = 0._r8
    rtmCTL%qsub(:,:)       = 0._r8
    rtmCTL%qgwl(:,:)       = 0._r8
    rtmCTL%discharge(:,:)  = 0._r8
    rtmCTL%volr(:)         = 0._r8
    rtmCTL%flood(:)        = 0._r8

  END SUBROUTINE RunoffInit

END MODULE RunoffMod
