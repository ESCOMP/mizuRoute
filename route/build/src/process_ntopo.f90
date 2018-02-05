module process_ntopo

! data types
USE nrtype                              ! variable types, etc.
USE nrtype,    only : integerMissing    ! missing value for integers
USE dataTypes, only : var_ilength       ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength       ! double precision type: var(:)%dat

! global vars
USE public_var                                ! public variables

! named variables
USE var_lookup,only:ixTOP,nVarsTOP            ! index of variables for the network topology

implicit none

! privacy -- everything private unless declared explicitly
private
public::ntopo

contains

 ! *********************************************************************
 ! public subroutine: read and process river network data
 ! *********************************************************************
 subroutine ntopo(&
                  ! output: model control
                  nHRU,             & ! number of HRUs
                  nSeg,             & ! number of stream segments
                  ! output: populate data structures
                  structHRU,        & ! ancillary data for HRUs
                  structSeg,        & ! ancillary data for stream segments
                  structHRU2seg,    & ! ancillary data for mapping hru2basin
                  structNTOPO,      & ! ancillary data for network toopology
                  ! output: error control
                  ierr, message)
 ! external subroutines : I/O
 use read_streamSeg,  only:getData               ! get the ancillary data
 use write_streamSeg, only:writeData             ! write the ancillary data
 ! external subroutines : network topology
 use network_topo,    only:hru2segment           ! get the mapping between HRUs and segments
 use network_topo,    only:up2downSegment        ! get the mapping between upstream and downstream segments
 use network_topo,    only:reach_list            ! reach list
 use network_topo,    only:reach_mask            ! identify all reaches upstream of a given reach
 ! This subroutine 1) read river network data and 2) populate river network topology data strucutres
 implicit none
 ! output: model control
 integer(i4b)      , intent(out)              :: nHRU             ! number of HRUs
 integer(i4b)      , intent(out)              :: nSeg             ! number of stream segments
 ! output: populate data structures
 type(var_dlength) , intent(out), allocatable :: structHRU(:)     ! HRU properties
 type(var_dlength) , intent(out), allocatable :: structSeg(:)     ! stream segment properties
 type(var_ilength) , intent(out), allocatable :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength) , intent(out), allocatable :: structNTOPO(:)   ! network topology
 ! output: error control
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)           :: cmessage           ! error message of downwind routine
 integer(i4b)                    :: tot_upstream       ! total number of all of the upstream stream segments for all stream segments
 integer(i4b)                    :: tot_upseg          ! total number of immediate upstream segments for all  stream segments
 integer(i4b)                    :: tot_hru            ! total number of all the upstream hrus for all stream segments
 integer(i4b)   , allocatable    :: ixHRU_desired(:)   ! indices of desired hrus
 integer(i4b)   , allocatable    :: ixSeg_desired(:)   ! indices of desired reaches
 integer*8                       :: time0,time1        ! times

 ! initialize error control
 ierr=0; message='ntopo/'

 ! initialize times
 call system_clock(time0)

 ! ---------- read in the stream segment information ---------------------------------------------------------

 ! get the number of HRUs and stream segments (needed for allocate statements)
 call getData(&
              ! input
              trim(ancil_dir)//trim(fname_ntop), & ! input: file name
              dname_nhru,   & ! input: dimension name of the HRUs
              dname_sseg,   & ! input: dimension name of the stream segments
              ! output: model control
              nHRU,         & ! output: number of HRUs
              nSeg,         & ! output: number of stream segments
              ! output: populate data structures
              structHRU,    & ! ancillary data for HRUs
              structSeg,    & ! ancillary data for stream segments
              structHRU2seg,& ! ancillary data for mapping hru2basin
              structNTOPO,  & ! ancillary data for network topology
              ! output: error control
              ierr,cmessage) ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after getData: time = ', time1-time0

 ! ---------- get the mapping between HRUs and segments ------------------------------------------------------

 ! get the mapping between HRUs and basins
 call hru2segment(&
                  ! input
                  nHRU,          & ! input: number of HRUs
                  nSeg,          & ! input: number of stream segments
                  ! input-output: data structures
                  structHRU,     & ! ancillary data for HRUs
                  structSeg,     & ! ancillary data for stream segments
                  structHRU2seg, & ! ancillary data for mapping hru2basin
                  structNTOPO,   & ! ancillary data for network toopology
                  ! output
                  tot_hru,    &   ! output: total number of all the upstream hrus for all stream segments
                  ierr, cmessage) ! output: error control

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after hru2segment: time = ', time1-time0

 ! ---------- get the mapping between upstream and downstream segments ---------------------------------------

 ! get the mapping between upstream and downstream segments
 call up2downSegment(&
                     ! input
                     nSeg,          & ! input: number of stream segments
                     ! input-output: data structures
                     structHRU,     & ! ancillary data for HRUs
                     structSeg,     & ! ancillary data for stream segments
                     structHRU2seg, & ! ancillary data for mapping hru2basin
                     structNTOPO,   & ! ancillary data for network toopology
                     ! output
                     tot_upseg,     & ! output: sum of immediate upstream segments
                     ierr, cmessage)  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after up2downSegment: time = ', time1-time0

 ! ---------- get the list of all upstream reaches above a given reach ---------------------------------------

 ! get the list of all upstream reaches above a given reach
 call reach_list(&
                 ! input
                 nSeg,                        & ! Number of reaches
                 structNTOPO,                 & ! Network topology
                 (computeReachList==compute), & ! flag to compute the reach list
                 ! output
                 tot_upstream,                & ! Total number of upstream reaches for all reaches
                 ierr, cmessage)                ! Error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after reach_list: time = ', time1-time0

 ! ---------- get indices of all segments above a prescribed reach ------------------------------------------

 ! identify all reaches upstream of a given reach
 call reach_mask(&
                 ! input
                 idSegOut,      &  ! input: reach index
                 structNTOPO,   &  ! input: network topology structures
                 nHRU,          &  ! input: number of HRUs
                 nSeg,          &  ! input: number of reaches
                 ! output: updated dimensions
                 tot_hru,       &  ! input+output: total number of all the upstream hrus for all stream segments
                 tot_upseg,     &  ! input+output: sum of immediate upstream segments
                 tot_upstream,  &  ! input+output: total number of upstream reaches for all reaches
                 ! output: dimension masks
                 ixHRU_desired, &  ! output: indices of desired hrus
                 ixSeg_desired, &  ! output: indices of desired reaches
                 ! output: error control
                 ierr, cmessage )  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 write(*,'(a,1x,i20)') 'after reach_mask: time = ', time1-time0

 print*, 'nDesire = ', size(ixHRU_desired)

 ! ---------- write network topology to a netcdf file -------------------------------------------------------

 ! write data
 call writeData(&
                ! input
                trim(ancil_dir)//trim(fname_output), & ! input: file name
                dname_nhru,    & ! input: dimension name of the HRUs
                dname_sseg,    & ! input: dimension name of the stream segments
                ! input: model control
                nHRU,          & ! input: number of HRUs
                nSeg,          & ! input: number of stream segments
                tot_hru,       & ! input: total number of all the upstream hrus for all stream segments
                tot_upseg,     & ! input: total number of immediate upstream segments for all  stream segments
                tot_upstream,  & ! input: total number of all of the upstream stream segments for all stream segments
                ! input: reach masks
                ixHRU_desired, & ! input: indices of desired hrus
                ixSeg_desired, & ! input: indices of desired reaches
                ! input: data structures
                structHRU,     & ! input: ancillary data for HRUs
                structSeg,     & ! input: ancillary data for stream segments
                structHRU2seg, & ! input: ancillary data for mapping hru2basin
                structNTOPO,   & ! input: ancillary data for network topology
                ! output: error control
                ierr,cmessage) ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif







 ! ---------- get the processing order -----------------------------------------------------------------------




 end subroutine ntopo

end module process_ntopo
