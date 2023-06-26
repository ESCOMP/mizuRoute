MODULE base_route

  ! Abstract class for routing method

  implicit none

  private
  public:: base_route_rch
  public:: routeContainer

  ! --- routing method container
  type :: routeContainer
    class(base_route_rch), allocatable :: rch_route
  end type


  ! --- routing method container
  type, abstract :: base_route_rch
    CONTAINS
    procedure(sub_route_rch), deferred :: route
  end type

  ABSTRACT INTERFACE

    SUBROUTINE sub_route_rch(this,          & !
                             iEns,          & ! input: ensemble index
                             segIndex,      & ! input: reach indice
                             ixDesire,      & ! input: index of verbose reach
                             T0, T1,        & ! input: start and end of simulation time step since start [sec]
                             NETOPO_in,     & ! input: reach topology data structure
                             RPARAM_in,     & ! input: reach parameter data structure
                             RCHSTA_out,    & ! inout: reach state data structure
                             RCHFLX_out,    & ! inout: reach flux data structure
                             ierr, message)   ! output: error control

      USE nrtype
      USE dataTypes, ONLY: STRFLX       ! fluxes in each reach
      USE dataTypes, ONLY: STRSTA       ! states in each reach
      USE dataTypes, ONLY: RCHTOPO      ! Network topology
      USE dataTypes, ONLY: RCHPRP       ! Reach parameter

      import base_route_rch
      ! Arguments
      class(base_route_rch)                     :: this
      integer(i4b),  intent(in)                 :: iEns             ! ensemble member
      integer(i4b),  intent(in)                 :: segIndex         ! ensemble member
      integer(i4b),  intent(in)                 :: ixDesire         ! index of the reach for verbose output
      real(dp),      intent(in)                 :: T0, T1           ! start and end of the time step (seconds)
      type(RCHTOPO), intent(in),    allocatable :: NETOPO_in(:)     ! River Network topology
      type(RCHPRP),  intent(inout), allocatable :: RPARAM_in(:)     ! River reach parameter
      type(STRSTA),  intent(inout)              :: RCHSTA_out(:,:)  ! reach state data
      type(STRFLX),  intent(inout)              :: RCHFLX_out(:,:)  ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
      integer(i4b),  intent(out)                :: ierr             ! error code
      character(*),  intent(out)                :: message          ! error message

    END SUBROUTINE sub_route_rch

  END INTERFACE

END MODULE base_route
