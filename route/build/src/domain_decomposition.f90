MODULE domain_decomposition

! External modules (general modules)
! numeric types
USE nrtype
! derived data types
USE dataTypes,         ONLY: var_ilength   ! integer type:     var(:)%dat
USE dataTypes,         ONLY: subbasin_omp  ! openMP domain data structure
! variable indices
USE var_lookup,        ONLY: ixNTOPO       ! index of variables for the netowork topolgy
! General utilities
USE nr_utility_module, ONLY: indexx        ! sorting array
USE nr_utility_module, ONLY: indexTrue     ! index at only true in array
USE nr_utility_module, ONLY: arth          ! generate sequential array
! updated and saved data
USE public_var

implicit none

private

public :: mpi_domain_decomposition
public :: omp_domain_decomposition

contains

 ! ***************************************************************
 ! public subroutine: MPI Domain decomposition
 ! ***************************************************************
 subroutine mpi_domain_decomposition(nNodes, nSeg, structNTOPO, nContribHRU, ierr, message)

   ! External modules
   USE globalData, only: domains                 ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain                 ! count of decomposed domains (tributaries + mainstems)

   implicit none
   ! Input variables
   integer(i4b),                   intent(in)  :: nNodes          ! number of nodes (root and computing nodes)
   integer(i4b),                   intent(in)  :: nSeg            ! number of stream segments
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)  ! network topology
   ! Output variables
   integer(i4b),                   intent(out) :: nContribHRU     ! total number of HRUs that are connected to a reach
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message         ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage        ! error message from subroutine

   ierr=0; message='mpi_domain_decomposition/'

   call classify_river_basin(nNodes,         &        ! input:  number of procs
                             nSeg,           &        ! input:  number of reaches in the entire river network
                             structNTOPO,    &        ! input:  river network data structure
                             domains,        &        ! output: domain data structure
                             nDomain,        &        ! output: number of domains
                             ierr, cmessage, &        ! output: error handling
                             nContribHRU=nContribHRU) ! output(optional): number of contributory HRUs for each reach
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call assign_node(nNodes, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine mpi_domain_decomposition

 ! ***************************************************************
 ! public subroutine: OMP domain decomposition
 ! ***************************************************************
 subroutine omp_domain_decomposition(nThreads, nSeg, structNTOPO, river_basin_out, ierr, message)

   ! External modules
   USE globalData, only: domains_omp             ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain_omp             ! count of decomposed domains (tributaries + mainstems)

   implicit none
   ! Input variables
   integer(i4b),                   intent(in)  :: nThreads        ! number of nodes (root and computing nodes)
   integer(i4b),                   intent(in)  :: nSeg            ! number of stream segments
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)  ! network topology
   ! Output variables
   type(subbasin_omp),allocatable, intent(out) :: river_basin_out(:)
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message         ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage        ! error message from subroutine

   ierr=0; message='omp_domain_decomposition/'

   call classify_river_basin(nThreads,       &        ! input:  number of threads
                             nSeg,           &        ! input:  number of reaches in the entire river network
                             structNTOPO,    &        ! input:  river network data structure
                             domains_omp,    &        ! output: domain data structure
                             nDomain_omp,    &        ! output: number of domains
                             ierr, cmessage)          ! output: error handling
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   call basin_order(nSeg, structNTOPO, river_basin_out, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine omp_domain_decomposition

 ! ***************************************************************
 ! private subroutine: Domain decomposition
 ! ***************************************************************
 subroutine classify_river_basin(nDivs,         & ! input:  number of divisions (nodes or threads)
                                 nSeg,          & ! input:  number of reaches in the entire river network
                                 structNTOPO,   & ! input:  river network data structure
                                 domains_out,   & ! output: domain data structure
                                 nDomains,      & ! output: number of domains
                                 ierr, message, & ! output: error handling
                                 nContribHRU)     ! output(optional): number of contributory HRUs for each reach
   ! Details:
   ! Divide the entire river basin, consiting of river reaches and HRUs, into tributaries (independent basins) and mainstems,
   ! such that the number of reaches in tributaries are less than threshold (= maxSegs). total upstream reaches above threshold
   ! are defined as mainstems.
   !
   ! The following data struct components need to be populated
   !   structNTOPO(:)%var(ixNTOPO%segIndex)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%downSegIndex)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%allUpSegIndices)%dat(:)
   !   structNTOPO(:)%var(ixNTOPO%nHRU)%dat(1)
   !   structNTOPO(:)%var(ixNTOPO%hruContribIx)%dat(:)
   !
   ! Populate the domain data structure
   !   domain(:)%basinType    : basin identifier 1 -> tributary, 2 -> mainstem
   !   domain(:)%segIndex(:)  : segment index within a basin
   !   domain(:)%hruIndex(:)  : hru indix within a basin
   !   domain(:)%idNode       : proc id (-1 through nNode-1) -1 is for mainstem but use pid=0

   ! updated and saved data
   USE public_var, only : maxDomain
   USE dataTypes,  only : subbasin_mpi  ! reach category (store mainstem code or pfaf code)

   implicit none

   ! Input variables
   integer(i4b),                   intent(in)  :: nDivs               ! number of nodes (root and computing nodes)
   integer(i4b),                   intent(in)  :: nSeg                ! number of stream segments
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)      ! network topology
   ! Output variables
   type(subbasin_mpi),             intent(out) :: domains_out(maxDomain)  ! domain decomposition data structure (maximum domain is set to maxDomain)
   integer(i4b),                   intent(out) :: nDomains
   integer(i4b),      optional,    intent(out) :: nContribHRU         ! total number of HRUs that are connected to a reach
   integer(i4b),                   intent(out) :: ierr
   character(len=strLen),          intent(out) :: message             ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage            ! error message from subroutine
   integer(i4b)                                :: maxSegs             ! upper limit of  number of tributary reaches
   integer(i4b),      allocatable              :: nHruLocal(:)        ! a number of HRU for selected reaches (e.g., reaches in one domain)
   integer(i4b)                                :: segIndex(nSeg)      ! reach index for all the reaches
   integer(i4b)                                :: downIndex(nSeg)     ! downstream reach index for all the reacheds
   logical(lgt)                                :: majorMainstem(nSeg) ! logical to indicate reach is "major" mainstem
   integer(i4b)                                :: nUpSeg              ! number of upstream reaches for a reach
   integer(i4b)                                :: sumHruLocal         ! sum of hrus that contribute to the segments
   integer(i4b)                                :: iSeg, ix            ! loop indices
   integer(i4b)                                :: ix1, ix2            ! first and last indices in array to subset
   logical(lgt)                                :: debug = .false.     ! print out reach info with node assignment for debugging

   integer*8                              :: cr, startTime, endTime
   real(dp)                               :: elapsedTime

   ierr=0; message='classify_river_basin/'
   call system_clock(count_rate=cr)

   ! check
   if (nSeg/=size(structNTOPO))then; ierr=20; message=trim(message)//'number of reach input is not correct'; return; endif

   ! Compute upper limit of reaches numbers within tributaries
   maxSegs = nSeg/nDivs

   ! Initialize number of domains
   nDomains = 0

   majorMainstem = .false.

   ! put ntopo data structure components into a separate array
   do iSeg = 1, nSeg
     segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
     downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
     nUpSeg          = size(structNTOPO(iSeg)%var(ixNTOPO%allUpSegIndices)%dat)
     if (nUpSeg > maxSegs) majorMainstem(iSeg) = .true.
   enddo

   call decomposeDomain(structNTOPO, majorMainstem, maxSegs, domains_out, nDomains, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

call system_clock(startTime)
   ! populate domain(:)%hruIndex
   if (present(nContribHRU)) nContribHRU = 0 ! total number of HRUs that contribute to the reach
   do ix=1,nDomains
     associate (ixSeg => domains_out(ix)%segIndex)
     allocate(nHruLocal(size(ixSeg)), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [nHruLocal]'; return; endif
     do iSeg = 1, size(ixSeg)
       nHruLocal(iSeg) = structNTOPO(ixSeg(iSeg))%var(ixNTOPO%nHRU)%dat(1)
     enddo
     sumHruLocal = sum(nHruLocal)
     if (present(nContribHRU)) nContribHRU=nContribHRU+sumHruLocal
     allocate(domains_out(ix)%hruIndex(sumHruLocal), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [domains_out(ix)%hruIndex]'; return; endif
     do iSeg = 1, size(ixSeg)
       sumHruLocal = sum(nHruLocal(1:iSeg))
       ix1 = sumHruLocal-nHruLocal(iSeg)+1
       ix2 = sumHruLocal
       domains_out(ix)%hruIndex(ix1:ix2) = structNTOPO(ixSeg(iSeg))%var(ixNTOPO%hruContribIx)%dat(:)
     enddo
     deallocate(nHruLocal, stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem deallocating [nHruLocal]'; return; endif
     end associate
   enddo
call system_clock(endTime)
elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
write(*,"(A,1PG15.7,A)") '   elapsed-time [hru_decomposition] = ', elapsedTime, ' s'

   if (debug) call print_screen()

   CONTAINS

   ! --------------------------------------------------
   !  FOR DEBUGGING
   ! --------------------------------------------------
   subroutine print_screen()
     ! debugging variables
     integer(i4b)                           :: segId(nSeg)            ! reach id for all the segments
     logical(lgt)                           :: missing(nSeg)

     do iSeg = 1,nSeg
       segId(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)
     end do

     print*,'seg_index segid down_index down_id node-id'
     do ix = 1,nDomains
      associate (segIndexSub => domains_out(ix)%segIndex, nSubSeg => size(domains_out(ix)%segIndex))
      do iSeg = 1,size(segIndexSub)
       if (downIndex(segIndexSub(iSeg)) > 0) then
       write(*,"(I9,A,I12,A,I9,A,I12,A,I2)") segIndexSub(iSeg),' ',segId(segIndexSub(iSeg)),' ', &
                                             downIndex(segIndexSub(iSeg)),' ',segId(downIndex(segIndexSub(iSeg))),' ', &
                                             domains_out(ix)%idNode
       else
       write(*,"(I9,A,I12,A,I9,A,I12,A,I2)") segIndexSub(iSeg),' ',segId(segIndexSub(iSeg)),' ', &
                                                     downIndex(segIndexSub(iSeg)),' ',-999,' ', &
                                                     domains_out(ix)%idNode
       endif
      end do
      end associate
     end do

     ! check not-assgined (missing) reaches
     missing = .true.
     do ix = 1,nDomains
      associate (segIndexSub => domains_out(ix)%segIndex)
      ! reach index array in order of node assignment
      do iSeg = 1,size(segIndexSub)
       missing(segIndexSub) = .false.
      end do
      end associate
     end do
     if (count(missing)>0) then
       print*,'segid down_id'
       do iSeg = 1, nSeg
        if (missing(iSeg)) then ! if not assigned reaches
         if (downIndex(iSeg)>0) then
          print*, segId(iSeg), segId(downIndex(iSeg))
         else
          print*, segId(iSeg), downIndex(iSeg)
         endif
        endif
       enddo
      else
       print*, 'NO MISSING SEGMENT: ALL SEGMENTS ARE ASSIGNED TO DOMAINS'
     endif

   end subroutine print_screen

 end subroutine classify_river_basin

 ! ***************************************************************
 ! private subroutine:
 ! ***************************************************************
 subroutine decomposeDomain(structNTOPO,  & ! input:
                            isMainstem,   & ! input:
                            maxSegs,      & ! input:
                            domains_out,  & ! inout: domain data structure
                            nDomains,     & ! inout: number of domains
                            ierr, message)
   ! Details:
   ! Given mainstem reaches, identify tributary basins (i.e., basin flows into mainstems)
   ! Update data structure -> domains
   ! domains(:)%basinType      code to indicate mainstem (0) or tributaries (1)
   !           %segIndex(:)    indices of reaches belong to this domain

   ! External modules
   USE dataTypes,          only: subbasin_mpi         ! reach category (store mainstem code or pfaf code)
   USE pfafstetter_module, only: lgc_tributary_outlet

   implicit none
   ! Input variables
   type(var_ilength), allocatable, intent(in)  :: structNTOPO(:)               ! network topology
   logical(lgt),                   intent(in)  :: isMainstem(:)                ! logical to indicate reach is mainstem
   integer(i4b),                   intent(in)  :: maxSegs                      ! threshold for upstream reach number to  define mainstem
   type(subbasin_mpi),             intent(inout) :: domains_out(maxDomain)     ! domain decomposition data structure (maximum domain is set to maxDomain)
   integer(i4b),                   intent(inout) :: nDomains                   ! number of domains (update)
   ! Output variables
   integer(i4b),                   intent(out) :: ierr                         ! error code
   character(len=strLen),          intent(out) :: message                      ! error message
   ! Local variables
   character(len=strLen)                       :: cmessage                     ! error message from subroutine
   integer(i4b)                                :: segIndex(size(structNTOPO))  ! reach indices for reaches within a basin
   integer(i4b)                                :: downIndex(size(structNTOPO)) ! downstream reach indices for reaches within a basin
   integer(i4b)                                :: iSeg,iTrib,iOut              ! loop indices
   integer(i4b)                                :: nSeg,nTrib,nOut,nMainstem    ! number of reaches, tributaries, basin outlets, and mainstems, respectively
   integer(i4b)                                :: nUpSegs                      ! numpber of upstream segments
   integer(i4b), allocatable                   :: ixOutlets(:)                 ! list of outlet reach indices for all the outlet reaches
   integer(i4b), allocatable                   :: ixTribOutlet(:)              ! tributary outlet reach indices
   integer(i4b), allocatable                   :: ixSubset(:)                  ! subset indices based on logical array from global index array
   logical(lgt), allocatable                   :: isTribOutlet(:)
   logical(lgt), allocatable                   :: isMainstem2d(:,:)
   integer(i4b), parameter                     :: mainstem=2
   integer(i4b), parameter                     :: tributary=1

   ierr=0; message='decomposeDomain/'

   ! put ntopo data structure components into a separate array
   nSeg = size(structNTOPO)
   do iSeg = 1, nSeg
     segIndex(iSeg)  = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)
     downIndex(iSeg) = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1)
   enddo

   ! Domain assignment for no-mainstem basin (i.e., basin size < maxSegs)
   ! Number of outlets
   nOut=count(downIndex<0)
   allocate(ixOutlets(nOut), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [ixOutlets]'; return; endif
   ixOutlets   = pack(segIndex, downIndex<0)
   do iOut = 1,nOut
     ! number of upstream reach for a tributary
     nUpSegs = size(structNTOPO(ixOutlets(iOut))%var(ixNTOPO%allUpSegIndices)%dat)
     if (nUpSegs > maxSegs) cycle
     nDomains = nDomains + 1
     domains_out(nDomains)%basinType = tributary
     domains_out(nDomains)%segIndex  = structNTOPO(ixOutlets(iOut))%var(ixNTOPO%allUpSegIndices)%dat
   enddo

   ! Count mainstem reaches and if there is no mainstem basin, exit.
   nMainstem = count(isMainstem)
   if (nMainstem == 0 ) return

   ! Domain assignment for mainstem basin (i.e., basin size > maxSegs) if exist
   ! 1. Populate mainstem segments
   ! Identify mainstem reach indices based on the basin reaches
   call indexTrue(isMainstem, ixSubset)

   ! 1. populate domains data structure
   nDomains = nDomains + 1
   allocate(domains_out(nDomains)%segIndex(nMainstem), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [domains_out(nDomains)%segIndex]'; return; endif
   domains_out(nDomains)%basinType = mainstem
   domains_out(nDomains)%segIndex = segIndex(ixSubset)

   deallocate(ixSubset, stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif

   ! 2. Populate tributary segments
   ! Get logical array indicating tributary outlet segments (size = nSeg)
   isMainstem2d = spread(isMainstem,2,1)  ! size: [nSeg x 2]
   call lgc_tributary_outlet(isMainstem2d, downIndex, isTribOutlet, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! Get indices of tributary outlet segments (size = number of tributaries)
   nTrib = count(isTribOutlet)

   ! Idenfity indices of tributary outlet reaches
   allocate(ixTribOutlet(nTrib), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [trib_outlet_idx]'; return; endif
   call indexTrue(isTribOutlet, ixSubset)
   ixTribOutlet = segIndex(ixSubset)
   deallocate(ixSubset, stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem deallocating [ixSubset]'; return; endif

   ! loop through each tributary to identify upstream reach indices
   do iTrib = 1,nTrib

     ! number of upstream reach for a tributary
     nUpSegs = size(structNTOPO(ixTribOutlet(iTrib))%var(ixNTOPO%allUpSegIndices)%dat)

     ! populate domains data structure
     ! increment number of domain
     nDomains = nDomains + 1

     allocate(domains_out(nDomains)%segIndex(nUpSegs), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [domains_out(nDomain)%segIndex]'; return; endif

     domains_out(nDomains)%basinType = tributary
     domains_out(nDomains)%segIndex = structNTOPO(ixTribOutlet(iTrib))%var(ixNTOPO%allUpSegIndices)%dat

   end do

 end subroutine decomposeDomain

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
   integer(i4b),              intent(in)  :: nNodes           ! nNodes
   ! Output variables
   integer(i4b),              intent(out) :: ierr
   character(len=strLen),     intent(out) :: message          ! error message
   ! Local variables
   character(len=strLen)                  :: cmessage         ! error message from subroutine
   integer(i4b)                           :: nTribSeg         ! number of tributary segments
   integer(i4b)                           :: nWork(nNodes-1) ! number of tributary workload for each node
   integer(i4b)                           :: ixNode(1)        ! node id where work load is minimum
   integer(i4b)                           :: ix,ixx           ! loop indices
   integer(i4b),allocatable               :: nSubSeg(:)       !
   integer(i4b),allocatable               :: rankDomain(:)   ! ranked domain based on size
   integer(i4b)                           :: nTrib            !
   integer(i4b)                           :: nEven            !
   integer(i4b)                           :: nSmallTrib       ! number of small tributaries to be processed in root node
   integer(i4b), parameter                :: tributary=1
   integer(i4b), parameter                :: mainstem=2
   logical(lgt),allocatable               :: isAssigned(:)

   ierr=0; message='assign_node/'

   allocate(nSubSeg(nDomain),rankDomain(nDomain),isAssigned(nDomain),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [nSubSeg,rankDomain,isAssigned]'; return; endif

   ! rank domain by number of segments
   ! count segments for each domain - nSubSeg
   ! count tributary domains - nTrib
   ! count total segments from tributary domains - nTribSeg
   ! rank domains based on nSubSeg - rankDomain

   nTribSeg = 0  ! number of tributary reaches
   nTrib = 0     ! number of tributaries
   do ix = 1,nDomain
    nSubSeg(ix) = size(domains(ix)%segIndex)
    if (domains(ix)%basinType==tributary) then
     nTribSeg = nTribSeg + nSubSeg(ix)
     nTrib = nTrib + 1
    endif
   end do
   call indexx(nSubSeg,rankDomain)

   ! root node is used to process small tributaries first in addition to mainstem later
   ! going through tributaries from the smallest, and accumulate number of tributary segments (nSmallTrib) up to "nEven"
   ! count a number of tributaries that are processed in mainstem cores

   nEven = nTribSeg/nNodes
   nSmallTrib=0
   isAssigned = .false.
   do ix = 1,nDomain
    ixx = rankDomain(ix)
    if (domains(ixx)%basinType==tributary) then
     nSmallTrib = nSmallTrib + nSubSeg(ixx)
     domains(ixx)%idNode = 0
     isAssigned(ixx) = .true.
     if(nSmallTrib > nEven) exit
    endif
   end do

   ! Distribute mainstem and "large" tributary to distributed cores (if more than one procs are available)
   nWork = 0
   do ix = nDomain,1,-1  ! Going through domain from the largest size
     ixx = rankDomain(ix)
     if (.not. isAssigned(ixx)) then
       if (domains(ixx)%basinType==mainstem) then   ! if domain is mainstem
         domains(ixx)%idNode = -1                   ! put -1 for temporarily but mainstem is handled in root proc (idNode = 0)
         isAssigned(ixx) = .true.
       elseif (domains(ixx)%basinType==tributary) then ! if domain is tributary
         ixNode = minloc(nWork)
         nWork(ixNode(1)) = nWork(ixNode(1))+size(domains(ixx)%segIndex)
         domains(ixx)%idNode = ixNode(1)
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
 ! private subroutine: Assign decomposed domain into procs
 ! ***************************************************************
 subroutine basin_order(nSeg, structNTOPO_in, river_basin_out, ierr, message)
   ! redo data structures

   ! External modules
   USE globalData, only: domains_omp       ! domain data structure - for each domain, pfaf codes and list of segment indices
   USE globalData, only: nDomain_omp       ! count of decomposed domains (tributaries + mainstems)

   implicit none
   ! Input variables
   integer(i4b),              intent(in)      :: nSeg
   type(var_ilength), allocatable, intent(in) :: structNTOPO_in(:)  ! network topology
   ! Output variables
   type(subbasin_omp),allocatable, intent(out) :: river_basin_out(:)!
   integer(i4b),              intent(out)     :: ierr
   character(len=strLen),     intent(out)     :: message            ! error message
   ! Local variables
   integer(i4b)                               :: segOrder(nSeg)     ! reach order
   integer(i4b)                               :: rankSegOrder(nSeg) ! ranked reach order
   integer(i4b)                               :: nTrib              ! number of tributary reaches
   integer(i4b)                               :: nMain              ! number of mainstem reaches
   integer(i4b)                               :: ix,ixx, iSeg       ! loop indices
   integer(i4b),allocatable                   :: nSubSeg(:)         !
   integer(i4b),allocatable                   :: subSegOrder(:)     ! reach order in a subset domain
   integer(i4b),allocatable                   :: rankDomain(:)      ! ranked domain based on size
   integer(i4b), parameter                    :: tributary=1
   integer(i4b), parameter                    :: mainstem=2
   logical(lgt),allocatable                   :: isAssigned(:)

   ierr=0; message='basin_order/'

   do iSeg=1,nSeg
     segOrder(iSeg) = structNTOPO_in(iSeg)%var(ixNTOPO%rchOrder)%dat(1)
   enddo

   ! sorting reach processing order
   call indexx(segOrder,rankSegOrder)

   allocate(river_basin_out(1), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [river_basin_out]'; return; endif

   allocate(nSubSeg(nDomain_omp),rankDomain(nDomain_omp),isAssigned(nDomain_omp),stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem allocating [nSubSeg,rankDomain,isAssigned]'; return; endif

   ! rank domains based on number of reaches i.e., nSubSeg - rankDomain
   ! count tributaries

   nTrib = 0; nMain = 0     ! initialize number of tributaries and mainstems
   do ix = 1,nDomain_omp
    nSubSeg(ix) = size(domains_omp(ix)%segIndex)
    if (domains_omp(ix)%basinType==tributary) nTrib=nTrib+1
    if (domains_omp(ix)%basinType==mainstem)  nMain=nMain+1
   end do
   call indexx(nSubSeg, rankDomain)

   if (nTrib/=0) then
     allocate(river_basin_out(1)%tributary(nTrib), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [river_basin_out%tributary]'; return; endif
   endif

   if (nMain/=0) then
     allocate(river_basin_out(1)%mainstem(nMain), stat=ierr)
     if(ierr/=0)then; message=trim(message)//'problem allocating [river_basin_out%mainstem]'; return; endif
   endif

   ! put reaches in tributaries and mainstem in the processing order within domain
   nTrib=0; nMain=0
   do ix = nDomain_omp,1,-1  ! Going through domain from the largest size

     ixx = rankDomain(ix)

     assigned:if (.not. isAssigned(ixx)) then

       associate(ixSegs => domains_omp(ixx)%segIndex)

       ! Compute reach order for only small basin
       allocate(subSegOrder(size(ixSegs)), stat=ierr)
       if(ierr/=0)then; message=trim(message)//'problem allocating segOrderTrib'; return; endif

       call indexx(rankSegOrder(ixSegs), subSegOrder)

       domain:if (domains_omp(ixx)%basinType==mainstem) then   ! if domain is mainstem

         nMain = nMain + 1

         allocate(river_basin_out(1)%mainstem(nMain)%segIndex(size(ixSegs)), stat=ierr)
         if(ierr/=0)then; message=trim(message)//'problem allocating river_basin_out(1)%mainstem(1)%segIndex'; return; endif

         river_basin_out(1)%mainstem(nMain)%segIndex(:) = ixSegs(subSegOrder)
         river_basin_out(1)%mainstem(nMain)%nRch        = size(ixSegs)

         isAssigned(ixx) = .true.

       elseif (domains_omp(ixx)%basinType==tributary) then ! if domain is tributary

         nTrib = nTrib + 1

         allocate(river_basin_out(1)%tributary(nTrib)%segIndex(size(ixSegs)), stat=ierr)
         if(ierr/=0)then; message=trim(message)//'problem allocating river_basin_out(1)%tributary(ix)%segIndex'; return; endif

         river_basin_out(1)%tributary(nTrib)%segIndex(:) = ixSegs(subSegOrder)
         river_basin_out(1)%tributary(nTrib)%nRch        = size(ixSegs)

         isAssigned(ixx) = .true.

       endif domain

       end associate

     endif assigned

   end do

 end subroutine basin_order

end module domain_decomposition
