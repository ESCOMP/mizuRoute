module process_ntopo

! data types
USE nrtype                              ! variable types, etc.
USE nrtype,    only : integerMissing    ! missing value for integers
USE dataTypes, only : var_ilength       ! integer type:          var(:)%dat
USE dataTypes, only : var_dlength       ! double precision type: var(:)%dat

! old data structures
USE nhru2basin, only : all_points             ! derived data types for hru2segment mapping
USE reachparam, only : RCHTOPO                ! Network topology structure
USE reachparam, only : RCHPRP                 ! Reach Parameter structure

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
                  ! output: data structures
                  hru2seg,          & ! hru-to-segment mapping structure
                  NETOPO,           & ! network topology structures
                  RPARAM,           & ! Reach Parameters
                  ! output: error control
                  ierr, message)
 ! external subroutines
 use read_streamSeg,    only:getData               ! get the ancillary data
 use read_streamSeg,    only:hru2segment           ! get the mapping between HRUs and segments
 use read_streamSeg,    only:up2downSegment        ! get the mapping between upstream and downstream segments
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
 ! output: data structures
 type(all_points)  , intent(out), allocatable :: hru2seg(:)       ! hru-segment mapping structure
 type(RCHTOPO)     , intent(out), allocatable :: NETOPO(:)        ! River Network topology
 type(RCHPRP)      , intent(out), allocatable :: RPARAM(:)        ! Reach Parameters
 ! output: error control
 integer(i4b)      , intent(out)              :: ierr             ! error code
 character(*)      , intent(out)              :: message          ! error message
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)          :: cmessage       ! error message of downwind routine
 integer(i4b)                   :: tot_upseg      ! total number of immediate upstream segments for all  stream segments
 integer(i4b)                   :: tot_hru        ! total number of all the upstream hrus for all stream segments
 integer*8                      :: time0,time1    ! times

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
              structNTOPO,  & ! ancillary data for network toopology
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
 write(*,'(a,1x,i20)') 'after hru2basin: time = ', time1-time0
     !
     !          ! ---------- get the mapping between upstream and downstream segments ---------------------------------------
     !
     !          ! get the mapping between upstream and downstream segments
     !          call up2downSegment(&
     !                              ! input
     !                              nSeg,         &    ! input: number of stream segments
     !                              sseg_acil,    &    ! input: stream segment parameters
     !                              ntop_acil,    &    ! input: network topology
     !                              hru2seg,      &    ! input: hru-segment mapping structure
     !                              ! output
     !                              NETOPO,       &    ! output: River Network topology
     !                              RPARAM,       &    ! output: Reach Parameters
     !                              tot_upseg,    &    ! output: sum of immediate upstream segments
     !                              ierr, cmessage)    ! output: error control
     !          if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     !
     !          ! get timing
     !          call system_clock(time1)
     !          write(*,'(a,1x,i20)') 'after assign_reachparam: time = ', time1-time0

 end subroutine ntopo

end module process_ntopo
