module remapping

  ! data types
  use nrtype
  use dataTypes, only : remap                  ! remapping data type
  use dataTypes, only : runoff                 ! runoff data type
  use dataTypes, only : var_ilength            ! integer type:          var(:)%dat
  use dataTypes, only : var_dlength            ! double precision type: var(:)%dat

  ! look-up variables
  use var_lookup,only:ixHRU,    nVarsHRU     ! index of variables for the HRUs
  use var_lookup,only:ixSEG,    nVarsSEG     ! index of variables for the stream segments
  use var_lookup,only:ixHRU2SEG,nVarsHRU2SEG ! index of variables for the hru2segment mapping
  use var_lookup,only:ixNTOPO,  nVarsNTOPO   ! index of variables for the network topology

  implicit none
  private
  public ::remap_runoff

  contains

  subroutine remap_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, message)
  implicit none
  ! input
  type(runoff)         , intent(in)  :: runoff_data      ! runoff for one time step for all HRUs
  type(remap)          , intent(in)  :: remap_data       ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(var_ilength)    , intent(in)  :: structHRU2seg(:) ! HRU-to-segment mapping
  ! output
  real(dp)             , intent(out) :: basinRunoff(:)   ! basin runoff
  integer(i4b)         , intent(out) :: ierr             ! error code
  character(len=strLen), intent(out) :: message          ! error message
  ! local
  integer(i4b)                       :: iHRU,jHRU        ! index of basin in the routing layer
  integer(i4b)                       :: ixOverlap        ! index in ragged array of overlapping polygons
  integer(i4b)                       :: ixRunoff         ! index in the runoff vector
  integer(i4b)                       :: ixPoly           ! loop through overlapping polygons for a given basin
  real(dp)                           :: sumWeights       ! used to check that the sum of weights equals one
  real(dp)    , parameter            :: xTol=1.e-6_dp    ! tolerance to avoid divide by zero
  integer(i4b), parameter            :: ixCheck=-huge(iHRU) ! basin to check
  !integer(i4b), parameter            :: ixCheck=24001479 ! basin to check

  ierr=0; message="remap_runoff/"

  ! initialize counter for the overlap vector
  ixOverlap = 1

  ! loop through basins in the routing layer
  do iHRU=1,size(structHRU2seg)

   ! define the HRU index in the routing vector
   jHRU = remap_data%hru_ix(iHRU)

   ! check that the basins match
   if( remap_data%hru_id(iHRU) /= structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1) )then
    message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
    ierr=20; return
   endif

   !print*, 'remap_data%hru_id(iHRU), structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1), remap_data%num_qhru(iHRU) = ', &
   !         remap_data%hru_id(iHRU), structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1), remap_data%num_qhru(iHRU)

   ! initialize the weighted average
   sumWeights        = 0._dp
   basinRunoff(jHRU) = 0._dp

   ! loop through the overlapping polygons
   do ixPoly=1,remap_data%num_qhru(iHRU) ! number of overlapping polygons

    ! check that the cell exists in the runoff file
    !print*, 'ixOverlap, remap_data%qhru_ix(ixOverlap) = ', ixOverlap, remap_data%qhru_ix(ixOverlap)
    if(remap_data%qhru_ix(ixOverlap)==integerMissing)then
     ixOverlap = ixOverlap + 1
     cycle
    endif

    ! get the index in the runoff file
    ixRunoff = remap_data%qhru_ix(ixOverlap)

    ! check that we have idenbtified the correct runoff HRU
    if( remap_data%qhru_id(ixOverlap) /= runoff_data%hru_id(ixRunoff) )then
     message=trim(message)//'mismatch in HRU ids for polygons in the runoff layer'
     ierr=20; return
    endif

    ! get the weighted average
    if(runoff_data%qSim(ixRunoff) > -xTol)then
     sumWeights        = sumWeights        + remap_data%weight(ixOverlap)
     basinRunoff(jHRU) = basinRunoff(jHRU) + remap_data%weight(ixOverlap)*runoff_data%qSim(ixRunoff)
    endif

    ! check
    if(remap_data%hru_id(iHRU)==ixCheck)then
     print*, 'remap_data%hru_id(iHRU)                         = ', remap_data%hru_id(iHRU)
     print*, 'structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1) = ', structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1)
     print*, 'remap_data%num_qhru(iHRU)                       = ', remap_data%num_qhru(iHRU)
     print*, 'ixRunoff, runoff_data%qSim(ixRunoff)            = ', ixRunoff, runoff_data%qSim(ixRunoff)
    endif

    !print*, 'remap_data%qhru_id(ixOverlap), runoff_data%hru_id(ixRunoff), remap_data%weight(ixOverlap), runoff_data%qSim(ixRunoff) = ', &
    !         remap_data%qhru_id(ixOverlap), runoff_data%hru_id(ixRunoff), remap_data%weight(ixOverlap), runoff_data%qSim(ixRunoff)

    ! increment the overlap index
    ixOverlap = ixOverlap + 1

   end do  ! looping through overlapping polygons

   ! compute weighted average
   if(sumWeights>xTol)then
    if(abs(1._dp - sumWeights)>xTol) basinRunoff(jHRU) = basinRunoff(jHRU) / sumWeights
   endif

   ! check
   if(remap_data%hru_id(iHRU)==ixCheck)then
    print*, 'basinRunoff(jHRU) = ', basinRunoff(jHRU)*86400._dp*1000._dp*365._dp
    print*, 'PAUSE : '; read(*,*)
   endif

   ! print progress
   if(mod(iHRU,100000)==0)then
    print*, trim(message)//'mapping runoff, iHRU, basinRunoff(jHRU) = ', &
                                            iHRU, basinRunoff(jHRU)
   endif

   !print*, 'basinRunoff(jHRU) = ', basinRunoff(jHRU)
   !print*, 'PAUSE : '; read(*,*)

  end do   ! looping through basins in the routing layer

  end subroutine remap_runoff
end module remapping
