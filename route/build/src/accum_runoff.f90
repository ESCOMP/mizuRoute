MODULE accum_runoff_module

! Accumulate upstream flow instantaneously
! This is not used as routed runoff.
! This is used to get total instantaneous upstream runoff at each reach

USE nrtype
USE dataTypes,   ONLY: STRFLX          ! fluxes in each reach
USE dataTypes,   ONLY: STRSTA          ! state in each reach
USE dataTypes,   ONLY: RCHTOPO         ! Network topology
USE dataTypes,   ONLY: RCHPRP          ! Reach parameter
USE public_var,  ONLY: iulog           ! i/o logical unit number
USE public_var,  ONLY: desireId        ! ID or reach where detailed reach state is print in log
USE globalData,  ONLY: idxSUM          ! routing method index for runoff accumulation method
USE base_route,  ONLY: base_route_rch  ! base (abstract) reach routing method class

implicit none

private
public::accum_runoff_rch

type, extends(base_route_rch) :: accum_runoff_rch
 CONTAINS
   procedure, pass :: route => accum_inst_runoff
end type accum_runoff_rch

CONTAINS

 ! *********************************************************************
 ! subroutine: perform accumulate immediate upstream flow
 ! *********************************************************************
 SUBROUTINE accum_inst_runoff(this,         &  ! "accum_runoff_rchr" object to bound this procedure
                              iEns,         &  ! input: index of runoff ensemble to be processed
                              segIndex,     &  ! input: index of reach to be processed
                              T0,T1,        &  ! input: start and end of the time step
                              NETOPO_in,    &  ! input: reach topology data structure
                              RPARAM_in,    &  ! input: reach parameter data structure
                              RCHSTA_out,   &  ! inout: reach state data structure
                              RCHFLX_out,   &  ! inout: reach flux data structure
                              ierr, message)   ! output: error control
 implicit none
 ! Argument variables
 class(accum_runoff_rch)                  :: this            ! "accum_runoff_rchr" object to bound this procedure
 integer(i4b),  intent(in)                :: iEns            ! runoff ensemble to be routed
 integer(i4b),  intent(in)                :: segIndex        ! segment where routing is performed
 real(dp),      intent(in)                :: T0,T1           ! start and end of the time step (seconds)
 type(RCHTOPO), intent(in),   allocatable :: NETOPO_in(:)    ! River Network topology
 type(RCHPRP),  intent(inout),allocatable :: RPARAM_in(:)    ! River reach parameter
 type(STRSTA),  intent(inout)             :: RCHSTA_out(:,:) ! reach state data
 type(STRFLX),  intent(inout)             :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 integer(i4b),  intent(out)               :: ierr            ! error code
 character(*),  intent(out)               :: message         ! error message
 ! Local variables
 real(dp)                                 :: q_upstream      ! upstream Reach fluxes
 integer(i4b)                             :: nUps            ! number of upstream segment
 integer(i4b)                             :: iUps            ! upstream reach index
 integer(i4b)                             :: iRch_ups        ! index of upstream reach in NETOPO
 character(len=strLen)                    :: fmt1,fmt2       ! format string

 ierr=0; message='accum_inst_runoff/'

 ! identify number of upstream segments of the reach being processed
 nUps = size(NETOPO_in(segIndex)%UREACHI)

 RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q = RCHFLX_out(iEns,segIndex)%BASIN_QR(1)

 q_upstream = 0._dp
 if (nUps>0) then

   do iUps = 1,nUps
     iRch_ups = NETOPO_in(segIndex)%UREACHI(iUps)      !  index of upstream of segIndex-th reach
     q_upstream = q_upstream + RCHFLX_out(iens,iRch_ups)%ROUTE(idxSUM)%REACH_Q
   end do

   RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q = RCHFLX_out(iEns,segIndex)%ROUTE(idxSUM)%REACH_Q + q_upstream

 endif

 ! check
 if(NETOPO_in(segIndex)%REACHID==desireId)then
   write(iulog,'(2a)') new_line('a'),'** Check upstream discharge accumulation **'
   write(iulog,'(a,1x,I10,1x,I10)') ' Reach index & ID =', segIndex, NETOPO_in(segIndex)%REACHID
   if (nUps>0) then
     write(fmt1,'(A,I5,A)') '(A,1X',nUps,'(1X,I10))'
     write(fmt2,'(A,I5,A)') '(A,1X',nUps,'(1X,F20.7))'
     write(iulog,'(a)')             ' * upstream reach index (NETOPO_in%UREACH) and discharge (uprflux) [m3/s] :'
     write(iulog,fmt1)              ' UREACHK =', (NETOPO_in(segIndex)%UREACHK(iUps), iUps=1,nUps)
     write(iulog,fmt2)              ' prflux  =', (RCHFLX_out(iens,NETOPO_in(segIndex)%UREACHI(iUps))%ROUTE(idxSUM)%REACH_Q, iUps=1,nUps)
   end if
   write(iulog,'(a)')             ' * local area discharge (RCHFLX_out%BASIN_QR(1)) and final discharge (RCHFLX_out%ROUTE(idxSUM)%REACH_Q) [m3/s] :'
   write(iulog,'(a,1x,G15.4)')     ' RCHFLX_out%BASIN_QR(1) =', RCHFLX_out(iEns,segIndex)%BASIN_QR(1)
   write(iulog,'(a,1x,G15.4)')     ' RCHFLX_out%ROUTE(idxSUM)%REACH_Q =', RCHFLX_out(iens,segIndex)%ROUTE(idxSUM)%REACH_Q
 endif

 END SUBROUTINE accum_inst_runoff

END MODULE accum_runoff_module
