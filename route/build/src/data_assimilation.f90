module data_assimilation

! general descriptions
!
!

USE nrtype
! data types
USE dataTypes, ONLY: STRFLX            ! fluxes in each reach
USE dataTypes, ONLY: STRSTA            ! state in each reach
USE dataTypes, ONLY: RCHTOPO           ! Network topology
! global data
USE public_var, ONLY: iulog             ! i/o logical unit number
USE public_var, ONLY: desireId          ! ID or reach where detailed reach state is print in log
USE public_var, ONLY: qBlendPeriod      ! number of time steps for which direct insertion is performed
USE public_var, ONLY: QerrTrend         ! temporal discharge error trend: 1->constant,2->linear, 3->logistic

implicit none

private
public::direct_insertion

contains

 ! *********************************************************************
 ! subroutine: perform direct-insertion
 ! *********************************************************************
 subroutine direct_insertion(iEns, segIndex,   & ! input: index of runoff ensemble to be processed
                             idxRoute,         & ! input: reachID to be checked by on-screen pringing
                             NETOPO_in,        & ! input: reach topology data structure
                             RCHSTA_out,       & ! inout: reach state data structure
                             RCHFLX_out,       & ! inout: reach flux data structure
                             ierr, message)      ! output: error control
   implicit none
   ! Argument variables
   integer(i4b),  intent(in)                 :: iEns              ! runoff ensemble to be routed
   integer(i4b),  intent(in)                 :: segIndex          ! segment where routing is performed
   integer(i4b),  intent(in)                 :: idxRoute          ! index of routing method
   type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)      ! River Network topology
   type(STRSTA),  intent(inout)              :: RCHSTA_out(:,:)   ! reach state data
   type(STRFLX),  intent(inout)              :: RCHFLX_out(:,:)   ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
   integer(i4b),  intent(out)                :: ierr              ! error code
   character(*),  intent(out)                :: message           ! error message
   ! Local variables
   logical(lgt)                              :: verbose           ! check details of variables
   real(dp)                                  :: Qcorrect          ! Discharge correction (when qmodOption=1) [m3/s]
   real(dp)                                  :: k                 ! the logistic decay rate or steepness of the curve.
   real(dp)                                  :: y0                !
   real(dp)                                  :: x0                !
   character(len=strLen)                     :: cmessage          ! error message from subroutine
   integer(i4b),parameter                    :: const=1           ! error reduction type: step function
   integer(i4b),parameter                    :: linear=2          ! error reduction type: linear function
   integer(i4b),parameter                    :: logistic=3        ! error reduction type: logistic function
   integer(i4b),parameter                    :: exponential=4     ! error reduction type: exponential function

   ierr=0; message='direct_insertion/'

   verbose = .false.
   if(NETOPO_in(segIndex)%REACHID == desireId) verbose = .true.

   if (RCHFLX_out(iens,segIndex)%Qobs>0._dp) then ! there is observation
     RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q - RCHFLX_out(iens,segIndex)%Qobs ! compute error
   end if

   if (RCHFLX_out(iens,segIndex)%Qelapsed > qBlendPeriod) then
     RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror=0._dp
   end if

   if (RCHFLX_out(iens,segIndex)%Qelapsed <= qBlendPeriod) then
     select case(QerrTrend)
       case(const)
         Qcorrect = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror
       case(linear)
         Qcorrect = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror*(1._dp - real(RCHFLX_out(iens,segIndex)%Qelapsed,dp)/real(qBlendPeriod, dp))
       case(logistic)
         x0 =0.25; y0 =0.90
         k = log(1._dp/y0-1._dp)/(qBlendPeriod/2._dp-qBlendPeriod*x0)
         Qcorrect = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror/(1._dp + exp(-k*(1._dp*RCHFLX_out(iens,segIndex)%Qelapsed-qBlendPeriod/2._dp)))
       case(exponential)
         if (RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror/=0._dp) then
           k = log(0.1_dp/abs(RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror))/(1._dp*qBlendPeriod)
           Qcorrect = RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%Qerror*exp(k*RCHFLX_out(iens,segIndex)%Qelapsed)
         else
           Qcorrect = 0._dp
         end if
       case default; message=trim(message)//'discharge error trend model must be 1(const),2(liear), or 3(logistic)'; ierr=81; return
     end select
   else
     Qcorrect=0._dp
   end if
   RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q = max(RCHFLX_out(iens,segIndex)%ROUTE(idxRoute)%REACH_Q-Qcorrect, 0._dp)

 end subroutine direct_insertion

end module data_assimilation
