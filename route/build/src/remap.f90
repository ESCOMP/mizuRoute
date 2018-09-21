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

  ! global data
  USE public_var,only:runoffMin, negRunoffTol, imiss
  USE globalData,only:time_conv,length_conv  ! conversion factors

  implicit none
  private
  public ::remap_runoff
  public ::basin2reach

  contains

  ! *****
  ! * public subroutine: used to map runoff data (on diferent grids/polygons) to the basins in the routing layer...
  ! ***************************************************************************************************************
  subroutine remap_runoff(runoff_data, remap_data, structHRU2seg, nSpace, basinRunoff, ierr, message)
  implicit none
  ! input
  type(runoff)         , intent(in)  :: runoff_data      ! runoff for one time step for all HRUs
  type(remap)          , intent(in)  :: remap_data       ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(var_ilength)    , intent(in)  :: structHRU2seg(:) ! HRU-to-segment mapping
  integer(i4b)         , intent(in)  :: nSpace(2)        ! size of spatial dimensions for runoff data
  ! output
  real(dp)             , intent(out) :: basinRunoff(:)   ! basin runoff
  integer(i4b)         , intent(out) :: ierr             ! error code
  character(len=strLen), intent(out) :: message          ! error message
  ! local
  character(len=strLen)              :: cmessage         ! error message from subroutine

  ierr=0; message="remap_runoff/"

  if (nSpace(2) == imiss) then
    call remap_1D_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
    call remap_2D_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  end subroutine remap_runoff

  ! *****
  ! private subroutine: used to map runoff data (on diferent polygons) to the basins in the routing layer...
  ! ***************************************************************************************************************
  subroutine remap_2D_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, message)
  implicit none
  ! input
  type(runoff)         , intent(in)  :: runoff_data         ! runoff for one time step for all HRUs
  type(remap)          , intent(in)  :: remap_data          ! data structure to remap data from a polygon (e.g., grid) to another polygon (e.g., basin)
  type(var_ilength)    , intent(in)  :: structHRU2seg(:)    ! HRU-to-segment mapping
  ! output
  real(dp)             , intent(out) :: basinRunoff(:)      ! basin runoff
  integer(i4b)         , intent(out) :: ierr                ! error code
  character(len=strLen), intent(out) :: message             ! error message
  ! local
  integer(i4b)                       :: iHRU,jHRU           ! index of basin in the routing layer
  integer(i4b)                       :: ixOverlap           ! index in ragged array of overlapping polygons
  integer(i4b)                       :: ii,jj               ! index of x and y in grid
  integer(i4b)                       :: ixPoly              ! loop through overlapping polygons for a given basin
  real(dp)                           :: sumWeights          ! used to check that the sum of weights equals one
  real(dp)    , parameter            :: xTol=1.e-6_dp       ! tolerance to avoid divide by zero
  integer(i4b), parameter            :: ixCheck=-huge(iHRU) ! basin to check
  integer(i4b), parameter            :: jxCheck=-huge(iHRU) ! basin to check
  logical(lgt), parameter            :: printWarn=.false.   ! flag to print warnings
  ierr=0; message="remap_2D_runoff/"

  ! initialize counter for the overlap vector
  ixOverlap = 1

  ! loop through hrus in the mapping layer
  do iHRU=1,size(remap_data%hru_ix)

   ! define the HRU index in the routing vector
   jHRU = remap_data%hru_ix(iHRU)

   ! if mapping data has hrus that do not exist in river network, skip that hru
   ! but increment index of weight and overlap-poly-id arrays
   if (jHRU == integerMissing)then
    if (remap_data%num_qhru(iHRU)/=integerMissing)then
      ixOverlap = ixOverlap + remap_data%num_qhru(iHRU)
    endif
    cycle
   endif

   ! check that the basins match
   if( remap_data%hru_id(iHRU) /= structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1) )then
    message=trim(message)//'mismatch in HRU ids for basins in the routing layer'
    ierr=20; return
   endif

   ! initialize the weighted average
   sumWeights        = 0._dp
   basinRunoff(jHRU) = 0._dp

   ! loop through the overlapping polygons
   do ixPoly=1,remap_data%num_qhru(iHRU) ! number of overlapping polygons

    ! index of i (x) and j (y) direction
    jj = remap_data%j_index(ixOverlap)
    ii = remap_data%i_index(ixOverlap)

    ! check i-indices
    if(ii < lbound(runoff_data%qSim2d,1) .or. ii > ubound(runoff_data%qSim2d,1))then
     if(printWarn) write(*,'(a,4(i0,a))') trim(message)//'WARNING: When computing weighted runoff at ', jHRU, 'th-HRU, i-index ', ii,' was not found in runoff grid data.'
     ixOverlap = ixOverlap + 1; cycle
    endif

    ! check j-indices
    if(jj < lbound(runoff_data%qSim2d,2) .or. jj > ubound(runoff_data%qSim2d,2))then
     if(printWarn) write(*,'(a,4(i0,a))') trim(message)//'WARNING: When computing weighted runoff at ', jHRU, 'th-HRU, j-index ', jj, 'was not found in runoff grid data.'
     ixOverlap = ixOverlap + 1; cycle
    endif

    ! get the weighted average
    if(runoff_data%qSim2d(ii,jj) > -xTol)then
     sumWeights        = sumWeights        + remap_data%weight(ixOverlap)
     basinRunoff(jHRU) = basinRunoff(jHRU) + remap_data%weight(ixOverlap)*runoff_data%qSim2D(ii,jj)
    endif

    ! check
    if(remap_data%i_index(iHRU)==ixCheck .and. remap_data%j_index(iHRU)==jxCheck)then
     print*, 'remap_data%i_index(iHRU),remap_data%j_index(iHRU) = ', remap_data%i_index(iHRU), remap_data%j_index(iHRU)
     print*, 'structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1)   = ', structHRU2seg(jHRU)%var(ixHRU2seg%hruId)%dat(1)
     print*, 'remap_data%num_qhru(iHRU)                         = ', remap_data%num_qhru(iHRU)
     print*, 'runoff_data%qSim2D(ii,jj)                         = ', runoff_data%qSim2D(ii,jj)
    endif

    ! increment the overlap index
    ixOverlap = ixOverlap + 1

   end do  ! looping through overlapping polygons

   ! compute weighted average
   if(sumWeights>xTol)then
    if(abs(1._dp - sumWeights)>xTol) basinRunoff(jHRU) = basinRunoff(jHRU) / sumWeights
   endif

   ! check
   if(remap_data%i_index(iHRU)==ixCheck .and. remap_data%j_index(iHRU)==jxCheck)then
    print*, 'basinRunoff(jHRU) = ', basinRunoff(jHRU)*86400._dp*1000._dp*365._dp
    print*, 'PAUSE : '; read(*,*)
   endif

   ! print progress
   !if(mod(iHRU,100000)==0)then
   ! print*, trim(message)//'mapping runoff, iHRU, basinRunoff(jHRU) = ', &
   !                                         iHRU, basinRunoff(jHRU)
   !endif

  end do   ! looping through basins in the mapping layer

  !!print*, 'PAUSE: after remap_2D_runoff'; read(*,*)

  end subroutine remap_2D_runoff

  ! *****
  ! private subroutine: used to map runoff data (on diferent polygons) to the basins in the routing layer...
  ! ***************************************************************************************************************
  subroutine remap_1D_runoff(runoff_data, remap_data, structHRU2seg, basinRunoff, ierr, message)
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

  ierr=0; message="remap_1D_runoff/"

  ! initialize counter for the overlap vector
  ixOverlap = 1

  ! loop through hrus in the mapping layer
  do iHRU=1,size(remap_data%hru_ix)

   ! define the HRU index in the routing vector
   jHRU = remap_data%hru_ix(iHRU)

   ! if mapping data has hrus that do not exist in river network, skip that hru
   ! but increment index of weight and overlap-poly-id arrays
   if (jHRU == integerMissing)then
    if (remap_data%num_qhru(iHRU)/=integerMissing)then
      ixOverlap = ixOverlap + remap_data%num_qhru(iHRU)
      print*, ixOverlap
    endif
    cycle
   endif

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
   !if(mod(iHRU,100000)==0)then
   ! print*, trim(message)//'mapping runoff, iHRU, basinRunoff(jHRU) = ', &
   !                                         iHRU, basinRunoff(jHRU)
   !endif

   !print*, 'basinRunoff(jHRU) = ', basinRunoff(jHRU)
   !print*, 'PAUSE : '; read(*,*)

  end do   ! looping through basins in the mapping layer

  end subroutine remap_1D_runoff

  ! *****
  ! * public subroutine: used to obtain streamflow for each stream segment...
  ! *************************************************************************

  subroutine basin2reach(&
                         ! input
                         basinRunoff,       & ! intent(in):  basin runoff (m/s)
                         structNTOPO,       & ! intent(in):  Network topology structure
                         structSEG,         & ! intent(in):  Network attributes structure
                         ! output
                         reachRunoff,       & ! intent(out): reach runoff (m/s)
                         ierr, message)       ! intent(out): error control
  implicit none
  ! input
  real(dp)             , intent(in)  :: basinRunoff(:)   ! basin runoff (m/s)
  type(var_ilength)    , intent(in)  :: structNTOPO(:)   ! Network topology structure
  type(var_dlength)    , intent(in)  :: structSEG(:)     ! Network attributes structure
  ! output
  real(dp)             , intent(out) :: reachRunoff(:)   ! reach runoff (m/s)
  integer(i4b)         , intent(out) :: ierr             ! error code
  character(len=strLen), intent(out) :: message          ! error message
  ! ----------------------------------------------------------------------------------------------
  ! local
  integer(i4b)                       :: iHRU             ! array index for contributing HRU
  integer(i4b)                       :: iSeg             ! array index for stream segment
  ! initialize error control
  ierr=0; message='basin2reach/'

  ! interpolate the data to the basins
  do iSeg=1,size(structSEG)

   ! associate variables in data structure
   associate(nContrib       => structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1),      & ! contributing HRUs
             hruContribIx   => structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat, & ! index of contributing HRU
             hruContribId   => structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat, & ! unique ids of contributing HRU
             basArea        => structSEG(  iSeg)%var(ixSEG%basArea)%dat(1),     & ! basin (total contributing HRU) area
             hruWeight      => structSEG(  iSeg)%var(ixSEG%weight)%dat          ) ! weight assigned to each HRU

   ! * case where HRUs drain into the segment
   if(nContrib > 0)then

    ! intialize the streamflow
    reachRunoff(iSeg) = 0._dp

    ! loop through the HRUs
    do iHRU=1,nContrib

     ! error check - runoff depth cannot be negative (no missing value)
     if( basinRunoff( hruContribIx(iHRU) ) < negRunoffTol )then
      write(message,'(a,i0)') trim(message)//'exceeded negative runoff tolerance for HRU ', hruContribId(iHRU)
      ierr=20; return
     endif

     ! compute the weighted average runoff depth (m/s)
     reachRunoff(iSeg) = reachRunoff(iSeg) + hruWeight(iHRU)*basinRunoff( hruContribIx(iHRU) )*time_conv*length_conv

    end do  ! (looping through contributing HRUs)

    ! ensure that routed streamflow is non-zero
    if(reachRunoff(iSeg) < runoffMin) reachRunoff(iSeg) = runoffMin

    ! convert basin average runoff volume (m3/s)
    reachRunoff(iSeg) = reachRunoff(iSeg)*basArea

   ! * special case where no HRUs drain into the segment
   else
    reachRunoff(iSeg) = runoffMin
   endif

   ! end association to data structures
   end associate

  end do  ! looping through stream segments

  end subroutine basin2reach

end module remapping
