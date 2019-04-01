MODULE domain_decomposition

! External modules (general modules)
! numeric types
USE nrtype
! derived data types
USE dataTypes,         ONLY: var_clength   ! character type:   var(:)%dat
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
! variable indices
USE var_lookup,        ONLY: ixPFAF        ! index of variables for the pfafstetter code
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
! General utilities
USE nr_utility_module, ONLY: indexx        ! sorting array
USE nr_utility_module, ONLY: indexTrue     ! I don't remember exactly
USE nr_utility_module, ONLY: arth          ! generate sequential array
! updated and saved data
USE public_var

implicit none

private

public :: classify_river_basin_omp
public :: classify_river_basin_mpi

contains

 ! ***************************************************************
 ! public subroutine: Main routine for MPI domain decomposition
 ! ***************************************************************
 subroutine classify_river_basin_mpi(nNodes,        & ! input:  number of procs
                                     nSeg,          & ! input:  number of reaches in the entire river network
                                     structPFAF,    & ! input:  pfafstetter code data structure
                                     structNTOPO,   & ! input:  river network data structure
                                     maxSegs,       & ! input:  threshold of tributary basin reaches
                                     nContribHRU,   & ! output: number of contributory HRUs for each reach
                                     ierr, message)   ! output: error handling
   ! Details:
   ! Divide the entire river basins, consiting of river reaches and HRUs, into tributaries (independent basins) and mainstems,
   ! such that the number of reaches in tributaries are less than threshold (= maxSegs). Using reach pfafstetter code to examine
   ! reaches/HRUs belong to tributaries/mainstems. Reaches missing pfafstetter code are grouped into one domain.
   !
   ! The following data struct components need to be populated
   !   structPFAF(:)%var(ixPFAF%code)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%segIndex)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%downSegIndex)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%allUpSegIndices)%dat(:)
   !   structNTOPO(:)%var(ixNTOPO%nHRU)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%hruContribIx)%dat(:)
   !
   ! Populate the domain data structure
   !   domain(:)%pfaf         : basin pfaf code
   !   domain(:)%segIndex(:)  : segment index within a basin
   !   domain(:)%hruIndex(:)  : hru indix within a basin
   !   domain(:)%isTrib       : T if a basin is a tributary, F if a basin is a mainstem
   !   domain(:)%idNode       : proc id (-1 through nNode-1) -1 is for mainstem but use pid=0

   ! updated and saved data
   USE globalData, only : domains                 ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only : nDomain                 ! count of decomposed domains (tributaries + mainstems)
   ! pfafstetter subroutines
   USE pfafstetter_module, only : get_common_pfaf

   implicit none

   ! Input variables
   integer(i4b),                   intent(in)  :: nNodes                 ! number of nodes (root and computing nodes)
   integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
   type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
   integer(i4b),                   intent(in)  :: maxSegs                ! threshold number of tributary reaches
   ! Output variables
   integer(i4b),                   intent(out) :: nContribHRU            ! total number of HRUs that are connected to a reach
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message                ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage               ! error message from subroutine
   character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the reaches
   character(len=32), allocatable              :: pfafs_sub(:)           ! pfaf_codes for a subset of the reaches
   character(len=32)                           :: pfafOutlet             ! pfaf_code for an outlet reach
   character(len=32)                           :: pfafCommon             ! common pfaf_codes over the entire basin
   character(len=32), allocatable              :: pfafOutlets(:)         ! list of pfaf_codes for all the outlet reaches
   integer(i4b)                                :: ixOutlet               ! reach index for an outlet reach
   integer(i4b),      allocatable              :: ixOutlets(:)           ! list of outlet reach indices for all the outlet reaches
   integer(i4b),      allocatable              :: ixSubset(:)            ! subset indices based on logical array from global index array
   integer(i4b),      allocatable              :: nHruLocal(:)           ! a number of HRU for selected reaches (e.g., reaches in one domain)
   integer(i4b)                                :: level                  ! number of digits of common pfaf codes given pfaf code at outlet reach
   integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the reaches
   integer(i4b)                                :: downIndex(nSeg)        ! downstream reach index for all the reacheds
   integer(i4b),      allocatable              :: segIndex_sub(:)        ! reach index for a subset of reacheds
   integer(i4b),      allocatable              :: downIndex_sub(:)       ! downstream reach index for a subset of reacheds
   integer(i4b)                                :: nOuts                  ! number of outlets
   integer(i4b)                                :: nInvalid               ! number of invalid pfafcode (0 or -999)
   integer(i4b)                                :: iSeg, iOut, ix         ! loop indices
   integer(i4b)                                :: ix1, ix2               ! first and last indices in array to subset
   logical(lgt),      allocatable              :: isInvalid(:)

   integer(i4b)                                :: segId(nSeg)            ! reach id for all the segments
   integer(i4b)                                :: ixSubSEG(nSeg)         !
   integer(i4b)                                :: ixSeg1, ixSeg2         !
   logical(lgt)                                :: seglgc(nSeg)

   ierr=0; message='classify_river_basin_mpi/'

   ! check
   if (nSeg/=size(structNTOPO))then; ierr=20; message=trim(message)//'number of reach input is not correct'; return; endif

   ! Initialize number of domains
   nDomain = 0

   ! put ntopo and pfaf data into a separate array
   forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
   forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
   forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)

   ! put segments with invalid pfaf code (0 or -999) into a separate domain
   nDomain = nDomain+1
   allocate(isInvalid(nSeg), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [isInvalid]'; return; endif
   isInvalid(1:nSeg)=.false.
   do iSeg = 1,nSeg
     if (trim(adjustl(pfafs(iSeg)))=='0' .or. trim(adjustl(pfafs(iSeg)))=='-999') then
       isInvalid(iSeg)=.true.
     end if
   end do
   nInvalid = count(isInvalid)
   allocate(ixSubset(nInvalid), domains(nDomain)%segIndex(nInvalid), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset,domains(nDomain)%segIndex]'; return; endif
   ixSubset = pack(arth(1,1,nSeg),isInvalid)
   domains(nDomain)%pfaf = '0'
   domains(nDomain)%segIndex = segIndex(ixSubset)

   ! Excluding invalid pfaf code reaches from the entire reaches
   allocate(segIndex_sub (nSeg-nInvalid), &
            downIndex_sub(nSeg-nInvalid), &
            pfafs_sub    (nSeg-nInvalid), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset]'; return; endif
   segIndex_sub  = pack(segIndex,(.not.isInvalid))
   downIndex_sub = pack(downIndex,(.not.isInvalid))
   pfafs_sub     = pack(pfafs,(.not.isInvalid))

   ! Number of outlets excluding invalid pfaf code basins
   nOuts=count(downIndex_sub<0)

   ! Outlet information - pfaf code and reach index
   allocate(pfafOutlets(nOuts), ixOutlets(nOuts), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [pfasOutlets, ixOutlet]'; return; endif
   pfafOutlets = pack(pfafs_sub, downIndex_sub<0)
   ixOutlets   = pack(segIndex_sub, downIndex_sub<0)

   ! process basin by basin
   do iOut = 1,nOuts

     ! pfaf code and reach index for this outlet
     pfafOutlet = adjustl(pfafOutlets(iOut))
     ixOutlet   = ixOutlets(iOut)

     ! get all the upstream reach indices for this outlet
     associate (ixUpSeg => structNTOPO(ixOutlet)%var(ixNTOPO%allUpSegIndices)%dat)

     ! Identify pfaf level for a river basin given outlet pfaf
     call get_common_pfaf(pfafs(ixUpSeg), pfafOutlet, level, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! perform domain decomposition for a river basin
     pfafCommon = trim(pfafOutlet(1:level))
     call decomposition(pfafs(ixUpSeg),     & ! input: basin-wide reach pfafcode
                        segIndex(ixUpSeg),  & ! input: basin-wide reach indices (index w.r.t. entire river network)
                        downIndex(ixUpSeg), & ! input: basin-wide downstream indices (index w.r.t. entire river network)
                        pfafCommon,         & ! input: common pfafcode for this basin (i.e., basin pfafcode)
                        maxSegs,            & ! input: max number of reaches for tributaries
                        ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end associate
   end do

   ! populate domain(:)%hruIndex
   nContribHRU = 0 ! total number of HRUs that contribute to the reach
   do ix=1,nDomain
     associate (ixSeg => domains(ix)%segIndex)
     allocate(nHruLocal(size(ixSeg)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [nHruLocal]'; return; endif
     do iSeg = 1, size(ixSeg)
       nHruLocal(iSeg) = structNTOPO(ixSeg(iSeg))%var(ixNTOPO%nHRU)%dat(1)
     enddo
     nContribHRU=nContribHRU+sum(nHruLocal)
     allocate(domains(ix)%hruIndex(sum(nHruLocal)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [domains(ix)%hruIndex]'; return; endif
     do iSeg = 1, size(ixSeg)
       ix1 = sum(nHruLocal(1:iSeg))-nHruLocal(iSeg)+1
       ix2 = sum(nHruLocal(1:iSeg))
       domains(ix)%hruIndex(ix1:ix2) = structNTOPO(ixSeg(iSeg))%var(ixNTOPO%hruContribIx)%dat(:)
     enddo
     deallocate(nHruLocal, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem deallocating [nHruLocal]'; return; endif
     end associate
   enddo

   ! Assign domains to node
   call assign_node(nNodes, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! --------------------------------------------------
   !  print to check
!   forall(iSeg=1:nSeg) segId(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)
!   print*,'segid global_index  pfaf-seg  pfaf-basin  basin-size  node-id'
!   ixSeg1=0; ixSeg2=0
!   ixSubSEG=-999
!   do ix = 1,nDomain
!    associate (segIndexSub => domains(ix)%segIndex, nSubSeg => size(domains(ix)%segIndex))
!    ! reach index array in order of node assignment
!    ixSeg1 = ixSeg2+1
!    ixSeg2 = ixSeg1+nSubSeg-1
!    ixSubSEG(ixSeg1:ixSeg2)  = domains(ix)%segIndex(1:nSubSeg)   ! global seg index per node
!    do iSeg = 1,size(segIndexSub)
!     write(*,"(I12,A,I9,A,A,A,A,A,I2)") segId(segIndexSub(iSeg)),' ',segIndexSub(iSeg),' ',trim(adjustl(pfafs(segIndexSub(iSeg)))),' ',trim(domains(ix)%pfaf),' ',domains(ix)%idNode
!    end do
!    end associate
!   end do
!   ! check not-assgined reaches
!   seglgc = .true.
!   do iSeg = 1, nSeg
!    if (ixSubSEG(iSeg)>0) then ! turned off at reach assigned properly
!     seglgc(ixSubSEG(iSeg)) = .false.
!    endif
!   enddo
!   do iSeg = 1, nSeg
!    if (seglgc(iSeg)) then ! check not assigned reaches
!     if (downIndex(iSeg)>0) then
!      print*, 'segid, downId, pfafcode = ', segId(iSeg), segId(downIndex(iSeg)), trim(pfafs(iSeg))
!     else
!      print*, 'segid, downId, pfafcode = ', segId(iSeg), downIndex(iSeg), trim(pfafs(iSeg))
!     endif
!    endif
!   enddo
!   ! --------------------------------------------------
!   print*,'----'
!   print*,'print out all the segments'
!   print*,'seq, segid, pfaf, down-segid, down-pfaf'
!   do ix = 1,nSeg
!     if (downIndex(ix)>0)then
!       write(*,"(I9,A,I9,A,A,A,I9,A,A)") ix,',',segId(ix),',',trim(adjustl(pfafs(ix))),',',segId(downIndex(ix)),',',trim(adjustl(pfafs(downIndex(ix))))
!     else
!       write(*,"(I9,A,I9,A,A,A,I9,A,A)") ix,',',segId(ix),',',trim(adjustl(pfafs(ix))),',',-999,',NaN'
!     endif
!   enddo

 end subroutine classify_river_basin_mpi


 ! ***************************************************************
 ! private subroutine: Assign decomposed domain into procs
 ! ***************************************************************
 subroutine assign_node(nNodes, ierr, message)
   ! assign domains into computing nodes

   ! External modules
   USE globalData, only: domains                 ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain                 ! count of decomposed domains (tributaries + mainstems)

   implicit none
   ! Input variables
   integer(i4b),                   intent(in)  :: nNodes           ! nNodes
   ! Output variables
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message          ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage         ! error message from subroutine
   integer(i4b)                                :: nTribSeg         ! number of tributary segments
   integer(i4b)                                :: nWork(nNodes-1)  ! number of tributary workload for each node
   integer(i4b)                                :: ixNode(1)        ! node id where work load is minimum
   integer(i4b)                                :: ix,ixx           ! loop indices
   integer(i4b),allocatable                    :: nSubSeg(:)       !
   integer(i4b),allocatable                    :: rnkRnDomain(:)   ! ranked domain based on size
   integer(i4b)                                :: nTrib            !
   integer(i4b)                                :: nEven            !
   integer(i4b)                                :: nSmallTrib       ! number of small tributaries to be processed in root node
   logical(lgt),allocatable                    :: isAssigned(:)

   ierr=0; message='assign_node/'

   allocate(nSubSeg(nDomain),rnkRnDomain(nDomain),isAssigned(nDomain),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [nSubSeg,rnkRnDomain,isAssigned]'; return; endif

   ! rank domain by number of segments
   ! count segments for each domain - nSubSeg
   ! count tributary domains - nTrib
   ! count total segments from tributary domains - nTribSeg
   ! rank domains based on nSubSeg - rnkRnDomain

   nTribSeg = 0  ! number of tributary reaches
   nTrib = 0     ! number of tributaries
   do ix = 1,nDomain
    nSubSeg(ix) = size(domains(ix)%segIndex)
    if (domains(ix)%pfaf(1:1)/='-') then
     nTribSeg = nTribSeg + nSubSeg(ix)
     nTrib = nTrib + 1
    endif
   end do
   call indexx(nSubSeg,rnkRnDomain)

   ! root node is used to process small tributaries first in addition to mainstem later
   ! going through tributaries from the smallest, and accumulate number of tributary segments (nSmallTrib) up to "nEven"
   ! count a number of tributaries that are processed in mainstem cores

   nEven = nTribSeg/nNodes
   nSmallTrib=0
   isAssigned = .false.
   do ix = 1,nDomain
    ixx = rnkRnDomain(ix)
    if (domains(ixx)%pfaf(1:1)/='-') then
     nSmallTrib = nSmallTrib + nSubSeg(ixx)
     domains(ixx)%idNode = 0
     domains(ixx)%isTrib = .true.
     isAssigned(ixx) = .true.
     if(nSmallTrib > nEven) exit
    endif
   end do

   ! Distribute mainstem and "large" tributary to distributed cores
   nWork(1:nNodes-1) = 0
   do ix = nDomain,1,-1  ! Going through domain from the largest size
    ixx = rnkRnDomain(ix)
    if (.not. isAssigned(ixx)) then
     if (domains(ixx)%pfaf(1:1)=='-') then   ! if domain is mainstem
      domains(ixx)%idNode = -1               ! put -1 for temporarily but mainstem is handled in root proc (idNode = 0)
      domains(ixx)%isTrib = .false.
      isAssigned(ixx) = .true.
     elseif (domains(ixx)%pfaf(1:1)/='-') then ! if domain is tributary
      ixNode = minloc(nWork)
      nWork(ixNode(1)) = nWork(ixNode(1))+size(domains(ixx)%segIndex)
      domains(ixx)%idNode = ixNode(1)
      domains(ixx)%isTrib = .true.
      isAssigned(ixx) = .true.
     endif
    endif
   end do

   ! check
   do ix = 1,nDomain
     if (.not.isAssigned(ix)) then
       write(cmessage, "(A,I1,A)") 'Domain ', ix, 'is not assigned to any nodes'
       ierr = 10; message=trim(message)//trim(cmessage); return
     endif
   enddo

 end subroutine assign_node


 ! ***************************************************************
 ! private subroutine: River network decomposition routine
 ! ***************************************************************
 recursive subroutine decomposition(pfafs, segIndex, downIndex, pfafCode, maxSegs, ierr, message)
   ! Details:
   ! Given a basin pfaf-code ("pfafCode") and reach information (reach pfaf-codes, reach index, and downstream index) for a particular basin,
   ! basin is further split into sub-basin domains with less than "maxSeg" of reaches
   !
   ! Given a basin pfaf-code (e.g.,968),
   ! 1.Examine number of reach, nSeg, for one higher level subbasin (e.g., 9681,9682,..9689)
   !
   ! 2.If the nSeg < maxSegs, call aggregation routine (see aggregation subroutine)
   !    2.1. if the basin is tributary or headwater, it is one domain, and update domains data structure
   !    2.2. if the basin is inter-basin, further this basin is decomposed into mainstems and tributaries and update domains data structure
   !
   ! 3. if the nSeg > maxSegs, call decomposition routine to examine one higher level subbasin
   !
   ! If a basin has less than threshold reaches, no decomposition process is performed and immediately populate data structure

   ! External modules
   USE globalData, only: domains                       ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain                       ! count of decomposed domains (tributaries + mainstems)

   implicit none

   ! Input variables
   character(len=32),              intent(in)  :: pfafs(:)        ! pfaf_code list
   integer(i4b),                   intent(in)  :: downIndex(:)    ! downstream segment index for all the segments
   integer(i4b),                   intent(in)  :: segIndex(:)     ! reach index for all the segments
   character(len=32),              intent(in)  :: pfafCode        ! common pfaf codes within a river basin
   integer(i4b)                                :: maxSegs         ! maximum number of reaches in basin
   ! Output variables
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message         ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage        ! error message from subroutine
   integer(i4b)                                :: iPfaf,iSeg,jSeg ! loop indices
   integer(i4b)                                :: nMatch          ! count for pfafs matching with pfafCode
   integer(i4b)                                :: nSeg            ! number of stream segments
   logical(lgt),     allocatable               :: isMatch(:)      ! logical to indicate xxxxx
   character(len=32)                           :: pfafCodeTmp     ! copy of pfafCode
   character(len=32),allocatable               :: subPfafs(:)     ! subset of pfaf_codes
   integer(i4b),     allocatable               :: subSegIndex(:)  ! subset of segment indices
   integer(i4b),     allocatable               :: subDownIndex(:) ! subset of donwstream segment original indices (indice are based on original segIndex)
   integer(i4b),     allocatable               :: downIndexNew(:) ! subset of downstream segment indices within subset of sgement indices
   character(len=32)                           :: subPfaf         ! pfaf_code appended by next level
   character(len=32)                           :: pfaf            ! a pfaf_code
   character(len=1)                            :: cPfaf           ! character digit

   ierr=0; message='decomposition/'

   nSeg = size(pfafs)

   ! if a basin has n reaches where n < maxSegs
   if (nSeg<=maxSegs) then
     nDomain = nDomain+1
     allocate(domains(nDomain)%segIndex(nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [domains(nDomain)%segIndex]'; return; endif
     ! populate domain data structures
     domains(nDomain)%pfaf = trim(pfafCode)
     domains(nDomain)%segIndex = segIndex
     return
   endif

   ! initialize sub-basin pfafcode
   pfafCodeTmp = trim(pfafCode)

   do iPfaf = 1,9
     write(cPfaf, '(I1)') iPfaf
     pfafCodeTmp = trim(pfafCodeTmp)//cPfaf

     if (.not. allocated(isMatch)) then
       allocate(isMatch(nSeg), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating [isMatch]'; return; endif
     endif
     isMatch(:) = .false.

     do iSeg = 1,size(pfafs)
       pfaf = adjustl(pfafs(iSeg))
       subPfaf = pfaf(1:len(trim(pfafCodeTmp)))
       if (subPfaf == pfafCodeTmp) then
         isMatch(iSeg) = .true.
       end if
     end do
     nMatch = count(isMatch)

     print*,'pfafCode, number of matches = ', trim(pfafCodeTmp), nMatch

     if (nMatch==0) then
       pfafCodeTmp = pfafCodeTmp(1:len(trim(pfafCodeTmp))-1)  ! decrement
       cycle
     end if

     if (nMatch < maxSegs) then
       print*,'Aggregate: nMatch less than nThreh = ', maxSegs
       allocate(subPfafs(nMatch), subDownIndex(nMatch), subSegIndex(nMatch), downIndexNew(nMatch), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating [subPfafs, subDownIndex, subSegIndex, downIndexNew]'; return; endif
       subPfafs     = pack(pfafs, isMatch)
       subSegIndex  = pack(segIndex, isMatch)
       subDownIndex = pack(downIndex, isMatch)
       ! redo donwstream index
       downIndexNew=-999
       do iSeg = 1, nMatch
         do jSeg = 1, nMatch
           if (subDownIndex(iSeg) == subSegIndex(jSeg)) then
             downIndexNew(iSeg) = jSeg
             exit
           end if
         end do
       end do
       call aggregate(pfafCodeTmp, subPfafs, subSegIndex, downIndexNew, ierr, cmessage) ! populate reach classification data structure
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
       deallocate(subPfafs, subDownIndex, subSegIndex, downIndexNew, stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem deallocating [subPfafs, subDownIndex, subSegIndex, downIndexNew]'; return; endif
     else
       print*,'Disaggregate: nMatch more than maxSegs = ', maxSegs
       call decomposition(pfafs, segIndex, downIndex, pfafCodeTmp, maxSegs, ierr, cmessage)
       if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     end if

     pfafCodeTmp = pfafCodeTmp(1:len(trim(pfafCodeTmp))-1) ! decrement

   end do

 end subroutine decomposition


 ! ***************************************************************
 ! private subroutine:
 ! ***************************************************************
 subroutine aggregate(pfafCode, pfafs, segIndex, downIndex, ierr, message)
   ! Details:
   ! Given a basin defined by "pfafCode",
   ! if a basin is tributaries, headwater, no further domain splits
   ! if a basin is inter-basin, split the basin into mainstem and tributaries.
   ! Assign minus pfafCode (pfafCode=968 -> -968) to mainstem reaches
   !
   ! Update data structure -> domains
   ! domains(:)%pfaf           code for mainstem reach or tributaries
   !           %segIndex(:)    indices of reaches that belong to pfaf

   ! External modules
   USE globalData, only: domains                       ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain                       ! count of decomposed domains (tributaries + mainstems)
   ! pfafstetter routines
   USE pfafstetter_module, only: lgc_tributary_outlet
   USE pfafstetter_module, only: find_mainstems
   USE pfafstetter_module, only: mainstem_code

   implicit none
   ! Input variables
   character(len=32),              intent(in)  :: pfafCode                  ! Common pfaf codes within a river basin
   character(len=32),              intent(in)  :: pfafs(:)                  ! pfaf codes within a basin
   integer(i4b),                   intent(in)  :: downIndex(:)              ! downstream reach indices for reaches within a basin
   integer(i4b),                   intent(in)  :: segIndex(:)               ! reach indices for reaches within a basin
   ! Output variables
   integer(i4b),                   intent(out) :: ierr                      ! error code
   character(len=strLen),          intent(out) :: message                   ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage                  ! error message from subroutine
   character(len=32)                           :: mainCode                  ! mainstem code
   character(len=32)                           :: subPfaf                   ! sub digit of segment pfaf code
   character(len=32)                           :: pfaf                      ! a segment pfafcode
   integer(i4b)                                :: iSeg,iTrib                ! loop indices
   integer(i4b)                                :: nSeg,nTrib,nMainstem      ! number of reaches, tributaries, and mainstems, respectively
   integer(i4b)                                :: nUpSegs                   ! numpber of upstream segments
   integer(i4b)                                :: pfafLen
   integer(i4b)                                :: pfaf_old, pfaf_new
   integer(i4b),allocatable                    :: ixSubset(:)               ! subset indices based on logical array from global index array
   character(len=32),allocatable               :: trib_outlet_pfafs(:)
   logical(lgt)                                :: isInterbasin, isHeadwater ! logical to indicate the basin is inter basin or headwater
   logical(lgt),allocatable                    :: isTrib(:)
   logical(lgt),allocatable                    :: isTribOutlet(:)
   logical(lgt),allocatable                    :: isMainstem(:)
   logical(lgt),allocatable                    :: isMainstem2d(:,:)

   ierr=0; message='aggregate/'

   ! Initialization
   nSeg = size(pfafs)

   ! get pfaf old and pfaf new
   pfafLen = len(trim(pfafCode))
   read(pfafCode(pfafLen-1:pfafLen-1),'(I1)') pfaf_old
   read(pfafCode(pfafLen:pfafLen),'(I1)') pfaf_new

   ! check if an interbasin or headwater
   isInterbasin = (mod(pfaf_new, 2)==1)
   isHeadwater  = (mod(pfaf_old, 2)==0 .and. pfaf_new==9)

   if (isInterbasin .and. (.not. isHeadwater)) then  ! if a river reach is in inter-basin and not headwater

     ! 1. Populate mainstem segments
     nDomain = nDomain + 1
     ! record the code for mainstem
     domains(nDomain)%pfaf = '-'//trim(pfafCode)
     ! get mainstem reaches (isMainstem(nSeg) logical array)
     call find_mainstems(pfafCode, pfafs, isMainstem, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     ! Count mainstem reaches and allocate for segIndex array
     nMainstem = count(isMainstem)
     allocate(domains(nDomain)%segIndex(nMainstem), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [domains(nDomain)%segIndex]'; return; endif
     ! Identify mainstem reach indices based on the basin reaches
     allocate(ixSubset(nMainstem), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset]'; return; endif
     ixSubset = pack(arth(1,1,nSeg),isMainstem)
     domains(nDomain)%segIndex = segIndex(ixSubset)
     deallocate(ixSubset, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif

     ! 2. Populate tributary segments
     ! Get logical array indicating tributary outlet segments (size = nSeg)
     isMainstem2d = spread(isMainstem,2,1)
     call lgc_tributary_outlet(isMainstem2d, downIndex, isTribOutlet, ierr, cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     ! Get indices of tributary outlet segments (size = number of tributaries)
     nTrib = count(isTribOutlet)
     ! idenfity tributary outlet reaches
     allocate(ixSubset(nTrib), trib_outlet_pfafs(nTrib), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset, trib_outlet_pfafs]'; return; endif
     ixSubset = pack(arth(1,1,nSeg),isTribOutlet)
     trib_outlet_pfafs = pfafs(ixSubset)
     deallocate(ixSubset, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif
     ! loop through each tributary to identify 1) tributary code (== mainstem code in tributary)
     !                                         2) reaches in tributary
     allocate(isTrib(nSeg), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [isTrib]'; return; endif
     do iTrib = 1,nTrib
       nDomain = nDomain + 1
       mainCode = mainstem_code(trim(trib_outlet_pfafs(iTrib)))
       domains(nDomain)%pfaf = mainCode

       isTrib(1:nSeg) = .false.
       do iSeg=1,nSeg
         pfaf = adjustl(pfafs(iSeg))
         subPfaf = pfaf(1:len(trim(mainCode)))
         if (trim(subPfaf) == trim(mainCode)) then
           isTrib(iSeg)=.true.
         endif
       enddo

       nUpSegs = count(isTrib)
       allocate(ixSubset(nUpSegs), domains(nDomain)%segIndex(nUpSegs), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating [ixSubset, domains(nDomain)%segIndex]'; return; endif
       ixSubset = pack(arth(1,1,nSeg),isTrib)
       domains(nDomain)%segIndex = segIndex(ixSubset)
       deallocate(ixSubset, stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif
     end do

   else   ! basin is tributaries, headwater

     nDomain = nDomain + 1
     allocate(domains(nDomain)%segIndex(size(segIndex)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating domains(nDomain)%segIndex for interbasin or headwater'; return; endif
     domains(nDomain)%pfaf = trim(pfafCode)
     domains(nDomain)%segIndex = segIndex

   end if

 end subroutine aggregate


 ! ***************************************************************
 ! public subroutine: Main routine for OMP domain decomposition
 ! ***************************************************************
 subroutine classify_river_basin_omp(nSeg, structPFAF, structNTOPO, river_basin, maxSegs, ierr, message)
  ! Identify tributary basin and mainstems using river network data and pfafstetter code
  ! Output: return populated basin dataType

  ! External modules
  ! derive data types
  USE dataTypes,          only: basin                ! basin data structure
  USE dataTypes,          only: reach                !
  ! pfafstetter routines
  USE pfafstetter_module, only: get_common_pfaf
  USE pfafstetter_module, only: lgc_mainstems
  USE pfafstetter_module, only: lgc_tributary_outlet
  USE pfafstetter_module, only: mainstem_code
  implicit none
  ! Input variables
  integer(i4b),                   intent(in)  :: nSeg                   ! number of stream segments
  type(var_clength), allocatable, intent(in)  :: structPFAF(:)          ! pfafstetter code
  type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)         ! network topology
  integer(i4b)                  , intent(in)  :: maxSegs                ! maximum reach number in each tributary
  ! Output variables
  type(basin),       allocatable, intent(out) :: river_basin(:)         ! river basin data
  integer(i4b),                   intent(out) :: ierr
  character(len=strLen),          intent(out) :: message                ! error message
  ! Local variables
  type(reach), allocatable                    :: tmpMainstem(:)         ! temporary reach structure for mainstem
  character(len=strLen)                       :: cmessage               ! error message from subroutine
  character(len=32)                           :: pfafs(nSeg)            ! pfaf_codes for all the segment
  character(len=32), allocatable              :: pfaf_out(:)            ! pfaf_codes for outlet segments
  character(len=32)                           :: upPfaf                 ! pfaf_code for one upstream segment
  character(len=32)                           :: outMainstemCode
  logical(lgt),      allocatable              :: mainstems(:,:)         ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: updated_mainstems(:,:) ! logical to indicate segment is mainstem at each level
  logical(lgt),      allocatable              :: lgc_trib_outlet(:)     ! logical to indicate segment is outlet of tributary
  logical(lgt)                                :: done                   ! logical
  integer(i4b)                                :: maxLevel=20
  integer(i4b)                                :: downIndex(nSeg)        ! downstream segment index for all the segments
  integer(i4b)                                :: segIndex(nSeg)         ! reach index for all the segments
  integer(i4b)                                :: segOrder(nSeg)         ! reach order for all the segments
  integer(i4b)                                :: rankSegOrder(nSeg)     ! ranked reach order for all the segments
  integer(i4b), allocatable                   :: segIndex_out(:)        ! index for outlet segment
  integer(i4b)                                :: level                  ! manstem level
  integer(i4b)                                :: nMains,nUpSegMain      ! number of mainstems in a level, and number of segments in a mainstem
  integer(i4b)                                :: nOuts                  ! number of outlets
  integer(i4b)                                :: nUpSegs                ! number of upstream segments for a specified segment
  integer(i4b)                                :: nTrib                  ! Number of Tributary basins
  integer(i4b), allocatable                   :: nSegTrib(:)            ! Number of segments in each tributary basin
  integer(i4b), allocatable                   :: nSegMain(:)            ! number of mainstem only segments excluding tributary segments
  integer(i4b), allocatable                   :: rankTrib(:)            ! index for ranked tributary based on num. of upstream reaches
  integer(i4b), allocatable                   :: rankMain(:)            ! index for ranked mainstem based on num. of mainstem reaches
  integer(i4b), allocatable                   :: segOrderTrib(:)        ! index for tributary upstream reach order
  integer(i4b), allocatable                   :: segMain(:)             ! index for segment in one mainstem
  integer(i4b), allocatable                   :: segOrderMain(:)        ! index for ordered segment in one mainstem
  integer(i4b), allocatable                   :: trPos(:)               ! tributary outlet indices
  integer(i4b), allocatable                   :: msPos(:)               ! mainstem indices
  integer(i4b)                                :: iSeg,jSeg,iOut         ! loop indices
  integer(i4b)                                :: iTrib,jTrib            ! loop indices
  integer(i4b)                                :: iMain,jMain            ! loop indices

  ierr=0; message='classify_river_basin_omp/'

  forall(iSeg=1:nSeg) pfafs(iSeg)     = structPFAF(iSeg)%var(ixPFAF%code)%dat(1)
  forall(iSeg=1:nSeg) downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
  forall(iSeg=1:nSeg) segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
  forall(iSeg=1:nSeg) segOrder(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)

  call indexx(segOrder,rankSegOrder)

  ! Number of outlets
  nOuts=count(downIndex<0)

  allocate(river_basin(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating river_basin'; return; endif
  allocate(pfaf_out(nOuts), segIndex_out(nOuts), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating [pfaf_out, segIndex_out]'; return; endif

  ! Outlet information - pfaf code and segment index
  pfaf_out = pack(pfafs, downIndex<0)
  segIndex_out = pack(segIndex, downIndex<0)

  do iOut = 1,nOuts

    print*, 'working on outlet:'//trim(adjustl(pfaf_out(iOut)))

    allocate(river_basin(iOut)%level(maxLevel), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem'; return; endif

    river_basin(iOut)%outIndex = segIndex_out(iOut)

!    call system_clock(startTime)
    ! Identify pfaf level given a river network
    call get_common_pfaf(pfafs, pfaf_out(iOut), level, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '       elapsed-time [get_common_pfaf] = ', elapsedTime, ' s'

!    call system_clock(startTime)
    ! Identify mainstem segments at all levels (up to maxLevel)
    call lgc_mainstems(pfafs, pfaf_out(iOut), maxLevel, mainstems, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_mainstems] = ', elapsedTime, ' s'

    ! Initial assignment of mainstem segments
    allocate(updated_mainstems(size(mainstems,1),size(mainstems,2)), stat=ierr)
    updated_mainstems = .false.
    updated_mainstems(:,level) = mainstems(:,level)

    ! Identify the lowest level mainstem segment index
    call indexTrue(mainstems(:,level), msPos)

    allocate(river_basin(iOut)%level(level)%mainstem(1), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%mainstem(level)%mainstem'; return; endif
    allocate(river_basin(iOut)%level(level)%mainstem(1)%segIndex(size(msPos)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%mainstem(:)%segIndex'; return; endif
    allocate(segOrderMain(size(msPos)), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem allocating segOrderMain'; return; endif

    ! Compute reach order for mainstem segment
    call indexx(rankSegOrder(msPos), segOrderMain)

    river_basin(iOut)%level(level)%mainstem(1)%segIndex(1:size(msPos))  = msPos(segOrderMain)
    river_basin(iOut)%level(level)%mainstem(1)%nRch = size(msPos)

!    call system_clock(startTime)
    ! Identify tributary outlets into a mainstem at the lowest level
    !  i.e. the segment that is not on any mainstems AND flows into any mainstem segments
    call lgc_tributary_outlet(mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
!    call system_clock(endTime)
!    elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!    write(*,"(A,1PG15.7,A)") '      elapsed-time [lgc_tributary_outlet] = ', elapsedTime, ' s'

    deallocate(mainstems, segIndex_out, segOrderMain, stat=ierr)
    if(ierr/=0)then; message=trim(message)//'problem deallocating mainstems'; return; endif

    do
      level = level + 1
      print*, 'Mainstem Level = ', level

      ! number of tributary basins
      nTrib = count(lgc_trib_outlet)

      allocate(nSegTrib(nTrib),rankTrib(nTrib), stat=ierr)
      if(ierr/=0)then; message=trim(message)//'problem allocating nSegTrib'; return; endif

      ! Extract array elements with only tributary outlet (keep indices in master array
      call indexTrue(lgc_trib_outlet, trPos)

      ! number of tributary segments
      do iTrib=1,nTrib
        nSegTrib(iTrib) = size(structNTOPO(trPos(iTrib))%var(ixNTOPO%allUpSegIndices)%dat)
      end do

      ! rank the tributary outlets based on number of upstream segments
      call indexx(nSegTrib, rankTrib)

      ! Identify mainstems of large tributaries (# of segments > maxSeg) and update updated_mainstems.
      done=.true.
      nMains = 0
      do iTrib=nTrib,1,-1
        if (nSegTrib(rankTrib(iTrib)) > maxSegs) then
           write(*,'(A,A,I5)') 'Exceed maximum number of segments: ', pfafs(trPos(rankTrib(iTrib))), nSegTrib(rankTrib(iTrib))
           nMains=nMains+1
           done=.false.
        else
           exit
        endif
      end do

      if (done) then ! if no mainstem/tributary updated, update tributary reach info, and then exist loop

!        print*, 'Main/Trib update finished, Now Update tributary reach ....'

        allocate(river_basin(iOut)%tributary(nTrib), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin%tributary'; return; endif

        do iTrib=nTrib,1,-1

!          call system_clock(startTime)

          jTrib = nTrib - iTrib + 1

!          if (mod(jTrib,1000)==0) print*, 'jTrib = ',jTrib

          allocate(river_basin(iOut)%tributary(jTrib)%segIndex(nSegTrib(rankTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(:)%tributary%segIndex'; return; endif
          allocate(segOrderTrib(nSegTrib(rankTrib(iTrib))), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating segOrderTrib'; return; endif

          ! compute reach order for tributary segments
          associate( segIndexTrib => structNTOPO(trPos(rankTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)
          call indexx(rankSegOrder(segIndexTrib), segOrderTrib)
          river_basin(iOut)%tributary(jTrib)%segIndex(:) = segIndexTrib(segOrderTrib)
          end associate

          ! compute number of segments in each tributary
          river_basin(iOut)%tributary(jTrib)%nRch = nSegTrib(rankTrib(iTrib))

          deallocate(segOrderTrib, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating segOrderTrib'; return; endif

!          if (mod(jTrib,1000)==0) then
!          call system_clock(endTime)
!          elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!          write(*,"(A,1PG15.7,A)") ' total elapsed-time [Update Trib. reach] = ', elapsedTime, ' s'
!          endif

        end do

        exit

      else ! if mainstem/tributary updated, store mainstem reach info at current level, update lgc_trib_outlet and go onto next level

        allocate(river_basin(iOut)%level(level)%mainstem(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%level(level)%mainstem'; return; endif
        allocate(tmpMainstem(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating tmpMainstem'; return; endif
        allocate(nSegMain(nMains), rankMain(nMains), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem allocating [nSegMain,rankMain]'; return; endif

        do iTrib=nTrib,nTrib-nMains+1,-1

!          if (mod(nTrib-iTrib+1,10)==0) print*, 'iTrib = ', nTrib-iTrib+1
!          call system_clock(startTime)

          ! Outlet mainstem code
          outMainstemCode = mainstem_code(pfafs(trPos(rankTrib(iTrib))))

          associate(upSegIndex => structNTOPO(trPos(rankTrib(iTrib)))%var(ixNTOPO%allUpSegIndices)%dat)

          nUpSegs = size(upSegIndex)

          allocate(segMain(nUpSegs), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating [segMain]'; return; endif

          nUpSegMain = 0
          do jSeg=1,nUpSegs
            upPfaf = pfafs(upSegIndex(jSeg))
            if (trim(mainstem_code(upPfaf)) /= trim(outMainstemCode)) cycle
            updated_mainstems(upSegIndex(jSeg),level) = .true.
            nUpSegMain = nUpSegMain+1
            segMain(nUpSegMain) = upSegIndex(jSeg)
          end do

          allocate(tmpMainstem(nTrib-iTrib+1)%segIndex(nUpSegMain), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating tmpMainstem(nUpSegMain)%segIndex'; return; endif
          allocate(segOrderMain(nUpSegMain), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating segOrderMain'; return; endif

          call indexx(rankSegOrder(segMain(1:nUpSegMain)), segOrderMain)
          tmpMainstem(nTrib-iTrib+1)%segIndex(:) = segMain(segOrderMain)
          nSegMain(nTrib-iTrib+1) = nUpSegMain

          end associate

          deallocate(segMain, segOrderMain, stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem deallocating [segMain, segOrderMain]'; return; endif

!          if (mod(nTrib-iTrib+1,10)==0) then
!          call system_clock(endTime)
!          elapsedTime = real(endTime-startTime, kind(dp))/10e8_dp
!          write(*,"(A,1PG15.7,A)") ' total elapsed-time [update trib each loop] = ', elapsedTime, ' s'
!          endif

        end do ! tributary loop

        ! rank the mainstem based on number of upstream segments (excluding tributaries)
        call indexx(nSegMain, rankMain)
        ! populate mainstem segment component in river basin structure
        do iMain=1,nMains
          jMain=rankMain(nMains-iMain+1)
          allocate(river_basin(iOut)%level(level)%mainstem(iMain)%segIndex(nSegMain(jMain)), stat=ierr)
          if(ierr/=0)then; message=trim(message)//'problem allocating river_basin(iOut)%level(level)%mainstem(iMain)%segIndex'; return; endif
          river_basin(iOut)%level(level)%mainstem(iMain)%segIndex(:) = tmpMainstem(jMain)%segIndex(:)
          river_basin(iOut)%level(level)%mainstem(iMain)%nRch = nSegMain(jMain)
        enddo

        ! update lgc_trib_outlet based on added mainstem
        call lgc_tributary_outlet(updated_mainstems(:,:level), downIndex, lgc_trib_outlet, ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

        deallocate(tmpMainstem, nSegMain, rankMain, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [tmpMainstem,nSegMain,rankMain]'; return; endif
        deallocate(nSegTrib, rankTrib, stat=ierr)
        if(ierr/=0)then; message=trim(message)//'problem deallocating [nSegTrib, rankTrib]'; return; endif

      endif
    end do ! end of tributary update loop

  end do ! outlet loop

 end subroutine classify_river_basin_omp


end module domain_decomposition
