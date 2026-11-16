! Lightweight pio_utils and globalData modules for isolating read_control and popMetadat in unit tests
MODULE pio_utils
  USE nrtype
  implicit none
  ! NetCDF type constants matching pio_utils.f90 for popMetadat initialization
  integer(i4b), parameter, public :: ncd_int    = 4 ! 32-bit integer NetCDF type code
  integer(i4b), parameter, public :: ncd_float  = 5 ! 32-bit float NetCDF type code
  integer(i4b), parameter, public :: ncd_double = 6 ! 64-bit double NetCDF type code
END MODULE pio_utils

MODULE globalData

  USE public_var
  USE nrtype
  USE dataTypes, ONLY: struct_info
  USE dataTypes, ONLY: dim_info
  USE dataTypes, ONLY: var_info
  USE objTypes,  ONLY: var_info_new
  USE var_lookup, ONLY: nStructures
  USE var_lookup, ONLY: nDimensions
  USE var_lookup, ONLY: nStateDims
  USE var_lookup, ONLY: nQdims
  USE var_lookup, ONLY: nVarsHRU
  USE var_lookup, ONLY: nVarsHRU2SEG
  USE var_lookup, ONLY: nVarsSEG
  USE var_lookup, ONLY: nVarsNTOPO
  USE var_lookup, ONLY: nVarsPFAF
  USE var_lookup, ONLY: nVarsRFLX
  USE var_lookup, ONLY: nVarsHFLX
  USE var_lookup, ONLY: nVarsBasinQ
  USE var_lookup, ONLY: nVarsIRFbas
  USE var_lookup, ONLY: nVarsBasTracer
  USE var_lookup, ONLY: nVarsIRF
  USE var_lookup, ONLY: nVarsKWT
  USE var_lookup, ONLY: nVarsKW
  USE var_lookup, ONLY: nVarsMC
  USE var_lookup, ONLY: nVarsDW
  USE var_lookup, ONLY: nVarsTracer

  implicit none

  save

  ! routing methods
  integer(i4b)                   , public :: nRoutes
  integer(i4b)                   , allocatable  , public :: routeMethods(:)
  logical(lgt)                   , public :: onRoute(0:nRouteMethods-1)
  integer(i4b)                   , public :: idxSUM
  integer(i4b)                   , public :: idxIRF
  integer(i4b)                   , public :: idxKWT
  integer(i4b)                   , public :: idxKW
  integer(i4b)                   , public :: idxMC
  integer(i4b)                   , public :: idxDW

  ! Misc. data
  character(len=strLen),           public :: runMode='standalone'
  logical(lgt),                    public :: isColdStart=.true.

  ! MPI/OMP variables
  logical(lgt),                    public :: masterproc=.false.

  ! conversion factors
  real(dp),                        public :: time_conv
  real(dp),                        public :: length_conv
  real(dp),                        public :: time_conv_solute
  real(dp),                        public :: mass_conv_solute

  ! general structure information
  type(struct_info),               public :: meta_struct(nStructures)
  type(dim_info),                  public :: meta_dims(nDimensions)
  type(dim_info),                  public :: meta_stateDims(nStateDims)
  type(dim_info),                  public :: meta_qDims(nQdims)
  type(dim_info),                  public :: meta_qDims_gage(nQdims)

  ! metadata structures
  type(var_info),                  public :: meta_HRU    (nVarsHRU    )
  type(var_info),                  public :: meta_HRU2SEG(nVarsHRU2SEG)
  type(var_info),                  public :: meta_SEG    (nVarsSEG    )
  type(var_info),                  public :: meta_NTOPO  (nVarsNTOPO  )
  type(var_info),                  public :: meta_PFAF   (nVarsPFAF   )
  type(var_info_new),              public :: meta_rflx   (nVarsRFLX   )
  type(var_info_new),              public :: meta_hflx   (nVarsHFLX   )
  type(var_info_new),              public :: meta_basinQ (nVarsBasinQ )
  type(var_info_new),              public :: meta_irf_bas(nVarsIRFbas )
  type(var_info_new),              public :: meta_bas_solute(nVarsBasTracer)
  type(var_info_new),              public :: meta_irf    (nVarsIRF    )
  type(var_info_new),              public :: meta_kwt    (nVarsKWT    )
  type(var_info_new),              public :: meta_kw     (nVarsKW     )
  type(var_info_new),              public :: meta_mc     (nVarsMC     )
  type(var_info_new),              public :: meta_dw     (nVarsDW     )
  type(var_info_new),              public :: meta_solute (nVarsTracer )

END MODULE globalData
