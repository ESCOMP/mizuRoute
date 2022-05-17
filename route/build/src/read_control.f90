MODULE read_control_module

USE nrtype
USE public_var

implicit none

private
public::read_control

CONTAINS

 ! =======================================================================================================
 ! * public subroutine: read the control file
 ! =======================================================================================================
 SUBROUTINE read_control(ctl_fname, err, message)

 ! global vars
 USE globalData, ONLY: time_conv,length_conv   ! conversion factors
 USE globalData, ONLY: masterproc              ! procs id and number of procs
 ! metadata structures
 USE globalData, ONLY: meta_HRU                ! HRU properties
 USE globalData, ONLY: meta_HRU2SEG            ! HRU-to-segment mapping
 USE globalData, ONLY: meta_SEG                ! stream segment properties
 USE globalData, ONLY: meta_NTOPO              ! network topology
 USE globalData, ONLY: meta_PFAF               ! pfafstetter code
 USE globalData, ONLY: meta_rflx               ! river flux variables
 USE globalData, ONLY: nRoutes                 ! number of active routing methods
 USE globalData, ONLY: routeMethods            ! active routing method index and id
 USE globalData, ONLY: onRoute                 ! logical to indicate actiive routing method(s)
 USE globalData, ONLY: idxSUM,idxIRF,idxKWT, &
                       idxKW,idxMC,idxDW
 ! named variables in each structure
 USE var_lookup, ONLY: ixHRU                   ! index of variables for data structure
 USE var_lookup, ONLY: ixHRU2SEG               ! index of variables for data structure
 USE var_lookup, ONLY: ixSEG                   ! index of variables for data structure
 USE var_lookup, ONLY: ixNTOPO                 ! index of variables for data structure
 USE var_lookup, ONLY: ixPFAF                  ! index of variables for data structure
 USE var_lookup, ONLY: ixRFLX, nVarsRFLX       ! index of variables for data structure
 ! external subroutines
 USE ascii_util_module, ONLY: file_open        ! open file (performs a few checks as well)
 USE ascii_util_module, ONLY: get_vlines       ! get a list of character strings from non-comment lines
 USE nr_utility_module, ONLY: char2int         ! convert integer number to a array containing individual digits

 implicit none
 ! argument variables
 character(*), intent(in)          :: ctl_fname               ! name of the control file
 integer(i4b),intent(out)          :: err                     ! error code
 character(*),intent(out)          :: message                 ! error message
 ! Local variables
 character(len=strLen),allocatable :: cLines(:)               ! vector of character strings
 character(len=strLen)             :: cName,cData             ! name and data from cLines(iLine)
 character(len=strLen)             :: cLength,cTime           ! length and time units
 integer(i4b)                      :: ipos                    ! index of character string
 integer(i4b)                      :: ibeg_name               ! start index of variable name in string cLines(iLine)
 integer(i4b)                      :: iend_name               ! end index of variable name in string cLines(iLine)
 integer(i4b)                      :: iend_data               ! end index of data in string cLines(iLine)
 integer(i4b)                      :: iLine                   ! index of line in cLines
 integer(i4b)                      :: iunit                   ! file unit
 integer(i4b)                      :: io_error                ! error in I/O
 integer(i4b)                      :: iRoute                  ! loop index
 character(len=strLen)             :: cmessage                ! error message from subroutine

 err=0; message='read_control/'

 ! *** get a list of character strings from non-comment lines ****
 ! open file (also returns un-used file unit used to open the file)
 call file_open(trim(ctl_fname),iunit,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

 ! get a list of character strings from non-comment lines
 call get_vlines(iunit,cLines,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

 close(iunit)

 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- read control file --- '
 end if

 ! loop through the non-comment lines in the input file, and extract the name and the information
 do iLine=1,size(cLines)

   ! initialize io_error
   io_error=0

   ! identify start and end of the name and the data
   ibeg_name = index(cLines(iLine),'<'); if(ibeg_name==0) err=20
   iend_name = index(cLines(iLine),'>'); if(iend_name==0) err=20
   iend_data = index(cLines(iLine),'!'); if(iend_data==0) err=20
   if(err/=0)then; message=trim(message)//'problem disentangling cLines(iLine) [string='//trim(cLines(iLine))//']';return;endif

   ! extract name of the information, and the information itself
   cName = adjustl(cLines(iLine)(ibeg_name:iend_name))
   cData = adjustl(cLines(iLine)(iend_name+1:iend_data-1))
   if (masterproc) then
     write(iulog,'(x,a,a,a)') trim(cName), ' --> ', trim(cData)
   endif

   ! populate variables
   select case(trim(cName))

   ! DIRECTORIES
   case('<ancil_dir>');            ancil_dir   = trim(cData)                           ! directory containing ancillary data
   case('<input_dir>');            input_dir   = trim(cData)                           ! directory containing input data
   case('<output_dir>');           output_dir  = trim(cData)                           ! directory containing output data
   case('<restart_dir>');          restart_dir = trim(cData)                           ! directory for restart output (netCDF)
   ! RUN CONTROL
   case('<case_name>');            case_name   = trim(cData)                           ! name of simulation. used as head of model output and restart file
   case('<sim_start>');            simStart    = trim(cData)                           ! date string defining the start of the simulation
   case('<sim_end>');              simEnd      = trim(cData)                           ! date string defining the end of the simulation
   case('<continue_run>');         read(cData,*,iostat=io_error) continue_run          ! logical; T-> append output in existing history files. F-> write output in new history file
   case('<newFileFrequency>');     newFileFrequency = trim(cData)                      ! frequency for new output files (day, month, annual, single)
   case('<route_opt>');            routOpt     = trim(cData)                           ! routing scheme options  0-> accumRunoff, 1->IRF, 2->KWT, 3-> KW, 4->MC, 5->DW
   case('<doesBasinRoute>');       read(cData,*,iostat=io_error) doesBasinRoute        ! basin routing options   0-> no, 1->IRF, otherwise error
   case('<seg_outlet>'   );        read(cData,*,iostat=io_error) idSegOut              ! desired outlet reach id (if -9999 --> route over the entire network)
   case('<is_lake_sim>');          read(cData,*,iostat=io_error) is_lake_sim           ! logical; lakes are simulated
   case('<is_flux_wm>');           read(cData,*,iostat=io_error) is_flux_wm            ! logical; provided fluxes to or from seg/lakes should be considered
   case('<is_vol_wm>');            read(cData,*,iostat=io_error) is_vol_wm             ! logical; provided target volume for managed lakes are considered
   case('<is_vol_wm_jumpstart>');  read(cData,*,iostat=io_error) is_vol_wm_jumpstart   ! logical; jump to the first time step target volume is set to true
   case('<suppress_runoff>');      read(cData,*,iostat=io_error) suppress_runoff       ! logical; suppress the read runoff to zero (0) no host model
   case('<suppress_P_Ep>');        read(cData,*,iostat=io_error) suppress_P_Ep         ! logical; suppress the precipitation/evaporation to zero (0) no host model
   ! RIVER NETWORK TOPOLOGY
   case('<fname_ntopOld>');        fname_ntopOld = trim(cData)                         ! name of file containing stream network topology information
   case('<ntopAugmentMode>');      read(cData,*,iostat=io_error) ntopAugmentMode       ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
   case('<fname_ntopNew>');        fname_ntopNew = trim(cData)                         ! name of file containing stream segment information
   case('<dname_nhru>');           dname_nhru    = trim(cData)                         ! dimension name of the HRUs
   case('<dname_sseg>');           dname_sseg    = trim(cData)                         ! dimension name of the stream segments
   ! RUNOFF, EVAPORATION AND PRECIPITATION FILE
   case('<fname_qsim>');           fname_qsim   = trim(cData)                          ! name of text file listing netcdf names. netCDF include runoff, evaporation and precipitation varialbes
   case('<vname_qsim>');           vname_qsim   = trim(cData)                          ! name of runoff variable
   case('<vname_evapo>');          vname_evapo  = trim(cData)                          ! name of actual evapoartion variable
   case('<vname_precip>');         vname_precip = trim(cData)                          ! name of precipitation variable
   case('<vname_time>');           vname_time   = trim(cData)                          ! name of time variable in the runoff file
   case('<vname_hruid>');          vname_hruid  = trim(cData)                          ! name of the HRU id
   case('<dname_time>');           dname_time   = trim(cData)                          ! name of time variable in the runoff file
   case('<dname_hruid>');          dname_hruid  = trim(cData)                          ! name of the HRU id dimension
   case('<dname_xlon>');           dname_xlon   = trim(cData)                          ! name of x (j,lon) dimension
   case('<dname_ylat>');           dname_ylat   = trim(cData)                          ! name of y (i,lat) dimension
   case('<units_qsim>');           units_qsim   = trim(cData)                          ! units of runoff
   case('<dt_qsim>');              read(cData,*,iostat=io_error) dt                    ! time interval of the gridded runoff
   case('<input_fillvalue>');      read(cData,*,iostat=io_error) input_fillvalue       ! fillvalue used for input variable
   ! FLUXES TO/FROM REACHES AND LAKE STATES FILE
   case('<fname_wm>');             fname_wm        = trim(cData)                       ! name of text file containing ordered nc file names
   case('<vname_flux_wm>');        vname_flux_wm   = trim(cData)                       ! name of varibale for fluxes to and from seg (reachs/lakes)
   case('<vname_vol_wm>');         vname_vol_wm    = trim(cData)                       ! name of varibale for target volume for managed lakes
   case('<vname_time_wm>');        vname_time_wm   = trim(cData)                       ! name of time variable
   case('<vname_segid_wm>');       vname_segid_wm  = trim(cData)                       ! name of the segid varibale in nc files
   case('<dname_time_wm>');        dname_time_wm   = trim(cData)                       ! name of time dimension
   case('<dname_segid_wm>');       dname_segid_wm  = trim(cData)                       ! name of the routing HRUs dimension
   ! RUNOFF REMAPPING
   case('<is_remap>');             read(cData,*,iostat=io_error) is_remap              ! logical case runnoff needs to be mapped to river network HRU
   case('<fname_remap>');          fname_remap          = trim(cData)                  ! name of runoff mapping netCDF
   case('<vname_hruid_in_remap>'); vname_hruid_in_remap = trim(cData)                  ! name of variable containing ID of river network HRU
   case('<vname_weight>');         vname_weight         = trim(cData)                  ! name of variable contating areal weights of runoff HRUs within each river network HRU
   case('<vname_qhruid>');         vname_qhruid         = trim(cData)                  ! name of variable containing ID of runoff HRU
   case('<vname_num_qhru>');       vname_num_qhru       = trim(cData)                  ! name of variable containing numbers of runoff HRUs within each river network HRU
   case('<vname_i_index>');        vname_i_index        = trim(cData)                  ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
   case('<vname_j_index>');        vname_j_index        = trim(cData)                  ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)
   case('<dname_hru_remap>');      dname_hru_remap      = trim(cData)                  ! name of dimension of river network HRU ID
   case('<dname_data_remap>');     dname_data_remap     = trim(cData)                  ! name of dimension of runoff HRU overlapping with river network HRU
   ! RESTART
   case('<restart_write>');        restart_write        = trim(cData)                  ! restart write option: N[n]ever, L[l]ast, S[s]pecified, Monthly, Daily
   case('<restart_date>');         restart_date         = trim(cData)                  ! specified restart date, yyyy-mm-dd (hh:mm:ss) for Specified option
   case('<restart_month>');        read(cData,*,iostat=io_error) restart_month         ! restart periodic month
   case('<restart_day>');          read(cData,*,iostat=io_error) restart_day           ! restart periodic day
   case('<restart_hour>');         read(cData,*,iostat=io_error) restart_hour          ! restart periodic hour
   case('<fname_state_in>');       fname_state_in       = trim(cData)                  ! filename for the channel states
   ! SPATIAL CONSTANT PARAMETERS
   case('<param_nml>');            param_nml            = trim(cData)                  ! name of namelist including routing parameter value
   ! USER OPTIONS: Define options to include/skip calculations
   case('<hydGeometryOption>');    read(cData,*,iostat=io_error) hydGeometryOption     ! option for hydraulic geometry calculations (0=read from file, 1=compute)
   case('<topoNetworkOption>');    read(cData,*,iostat=io_error) topoNetworkOption     ! option for network topology calculations (0=read from file, 1=compute)
   case('<computeReachList>');     read(cData,*,iostat=io_error) computeReachList      ! option to compute list of upstream reaches (0=do not compute, 1=compute)
   ! TIME
   case('<time_units>');           time_units = trim(cData)                            ! time units. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
   case('<calendar>');             calendar   = trim(cData)                            ! calendar name
   ! GAUGE META
   case('<gageMetaFile>');         gageMetaFile = trim(cData)                          ! name of csv file containing gauge metadata (gauge id, reach id, gauge lat/lon)
   case('<gageOnlyOutput>');       read(cData,*,iostat=io_error) gageOnlyOutput        ! logical; T-> history file output at only gauge points
   ! MISCELLANEOUS
   case('<debug>');                read(cData,*,iostat=io_error) debug                 ! print out detailed information throught the probram
   case('<desireId>'   );          read(cData,*,iostat=io_error) desireId              ! turn off checks or speficy reach ID if necessary to print on screen
   ! PFAFCODE
   case('<maxPfafLen>');           read(cData,*,iostat=io_error) maxPfafLen            ! maximum digit of pfafstetter code (default 32)
   case('<pfafMissing>');          pfafMissing = trim(cData)                           ! missing pfafcode (e.g., reach without any upstream area)
   ! OUTPUT OPTIONS
   case('<basRunoff>');            read(cData,*,iostat=io_error) meta_rflx(ixRFLX%basRunoff        )%varFile  ! default: true
   case('<instRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%instRunoff       )%varFile  ! default: false
   case('<dlayRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%dlayRunoff       )%varFile  ! default: false
   case('<sumUpstreamRunoff>');    read(cData,*,iostat=io_error) meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile  ! default: false
   case('<KWTroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWTroutedRunoff  )%varFile  ! default: true (turned off if inactive)
   case('<IRFroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%IRFroutedRunoff  )%varFile  ! default: true (turned off if inactive)
   case('<KWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<DWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<MCroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCroutedRunoff   )%varFile  ! default: true (turned off if inactive)
   case('<volume>');               read(cData,*,iostat=io_error) meta_rflx(ixRFLX%volume           )%varFile  ! default: true

   ! VARIABLE NAMES for data (overwrite default name in popMeta.f90)
   ! HRU structure
   case('<varname_area>'         ); meta_HRU    (ixHRU%area            )%varName = trim(cData) ! HRU area

   ! Mapping from HRUs to stream segments
   case('<varname_HRUid>'        ); meta_HRU2SEG(ixHRU2SEG%HRUid       )%varName = trim(cData) ! HRU id
   case('<varname_HRUindex>'     ); meta_HRU2SEG(ixHRU2SEG%HRUindex    )%varName = trim(cData) ! HRU index
   case('<varname_hruSegId>'     ); meta_HRU2SEG(ixHRU2SEG%hruSegId    )%varName = trim(cData) ! the stream segment id below each HRU
   case('<varname_hruSegIndex>'  ); meta_HRU2SEG(ixHRU2SEG%hruSegIndex )%varName = trim(cData) ! the stream segment index below each HRU

   ! reach properties
   case('<varname_length>'         ); meta_SEG    (ixSEG%length          )%varName =trim(cData)   ! length of segment  (m)
   case('<varname_slope>'          ); meta_SEG    (ixSEG%slope           )%varName =trim(cData)   ! slope of segment   (-)
   case('<varname_width>'          ); meta_SEG    (ixSEG%width           )%varName =trim(cData)   ! width of segment   (m)
   case('<varname_man_n>'          ); meta_SEG    (ixSEG%man_n           )%varName =trim(cData)   ! Manning's n        (weird units)
   case('<varname_hruArea>'        ); meta_SEG    (ixSEG%hruArea         )%varName =trim(cData)   ! local basin area (m2)
   case('<varname_weight>'         ); meta_SEG    (ixSEG%weight          )%varName =trim(cData)   ! HRU weight
   case('<varname_timeDelayHist>'  ); meta_SEG    (ixSEG%timeDelayHist   )%varName =trim(cData)   ! time delay histogram for each reach (s)
   case('<varname_upsArea>'        ); meta_SEG    (ixSEG%upsArea         )%varName =trim(cData)   ! area above the top of the reach -- zero if headwater (m2)
   case('<varname_basUnderLake>'   ); meta_SEG    (ixSEG%basUnderLake    )%varName =trim(cData)   ! Area of basin under lake  (m2)
   case('<varname_rchUnderLake>'   ); meta_SEG    (ixSEG%rchUnderLake    )%varName =trim(cData)   ! Length of reach under lake (m)
   case('<varname_minFlow>'        ); meta_SEG    (ixSEG%minFlow         )%varName =trim(cData)   ! minimum environmental flow
   case('<varname_D03_MaxStorage>' ); meta_SEG    (ixSEG%D03_MaxStorage  )%varName =trim(cData)   ! Doll 2006; maximum active storage for Doll 2003 formulation
   case('<varname_D03_Coefficient>'); meta_SEG    (ixSEG%D03_Coefficient )%varName =trim(cData)   ! Doll 2006; coefficient for Doll 2003 formulation (day-1)
   case('<varname_D03_Power>'      ); meta_SEG    (ixSEG%D03_Power       )%varName =trim(cData)   ! Doll 2006; power for Doll 2003 formulation

   case('<varname_HYP_E_emr>'      ); meta_SEG    (ixSEG%HYP_E_emr       )%varName =trim(cData)   ! HYPE; elevation of emergency spillway [m]
   case('<varname_HYP_E_lim>'      ); meta_SEG    (ixSEG%HYP_E_lim       )%varName =trim(cData)   ! HYPE; elevation below which primary spillway flow is restrcited [m]
   case('<varname_HYP_E_min>'      ); meta_SEG    (ixSEG%HYP_E_min       )%varName =trim(cData)   ! HYPE; elevation below which outflow is zero [m]
   case('<varname_HYP_E_zero>'     ); meta_SEG    (ixSEG%HYP_E_zero      )%varName =trim(cData)   ! HYPE; elevation at which lake/reservoir storage is zero [m]
   case('<varname_HYP_Qrate_emr>'  ); meta_SEG    (ixSEG%HYP_Qrate_emr   )%varName =trim(cData)   ! HYPE; emergency rate of flow for each unit of elevation above HYP_E_emr [m3/s]
   case('<varname_HYP_Erate_emr>'  ); meta_SEG    (ixSEG%HYP_Erate_emr   )%varName =trim(cData)   ! HYPE; power for the rate of flow for each unit of elevation above HYP_E_emr [-]
   case('<varname_HYP_Qrate_prim>' ); meta_SEG    (ixSEG%HYP_Qrate_prim  )%varName =trim(cData)   ! HYPE; the average yearly or long term output from primary spillway [m3/s]
   case('<varname_HYP_Qrate_amp>'  ); meta_SEG    (ixSEG%HYP_Qrate_amp   )%varName =trim(cData)   ! HYPE; amplitude of the Qrate_main [-]
   case('<varname_HYP_Qrate_phs>'  ); meta_SEG    (ixSEG%HYP_Qrate_phs   )%varName =trim(cData)   ! HYPE; phase of the Qrate_main based on the day of the year [-]; default 100
   case('<varname_HYP_prim_F>'     ); meta_SEG    (ixSEG%HYP_prim_F      )%varName =trim(cData)   ! HYPE; if the reservoir has a primary spillway then set to 1 otherwise 0
   case('<varname_HYP_A_avg>'      ); meta_SEG    (ixSEG%HYP_A_avg       )%varName =trim(cData)   ! HYPE; average area for the lake; this might not be used if bathymetry is provided [m]

   case('<varname_H06_Smax>'       ); meta_SEG    (ixSEG%H06_Smax        )%varName =trim(cData)   ! Hanasaki 2006; maximume reservoir storage [m3]
   case('<varname_H06_alpha>'      ); meta_SEG    (ixSEG%H06_alpha       )%varName =trim(cData)   ! Hanasaki 2006; fraction of active storage compared to total storage [-]
   case('<varname_H06_envfact>'    ); meta_SEG    (ixSEG%H06_envfact     )%varName =trim(cData)   ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
   case('<varname_H06_S_ini>'      ); meta_SEG    (ixSEG%H06_S_ini       )%varName =trim(cData)   ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
   case('<varname_H06_c1>'         ); meta_SEG    (ixSEG%H06_c1          )%varName =trim(cData)   ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
   case('<varname_H06_c2>'         ); meta_SEG    (ixSEG%H06_c2          )%varName =trim(cData)   ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
   case('<varname_H06_exponent>'   ); meta_SEG    (ixSEG%H06_exponent    )%varName =trim(cData)   ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
   case('<varname_H06_denominator>'); meta_SEG    (ixSEG%H06_denominator )%varName =trim(cData)   ! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
   case('<varname_H06_c_compare>'  ); meta_SEG    (ixSEG%H06_c_compare   )%varName =trim(cData)   ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
   case('<varname_H06_frac_Sdead>' ); meta_SEG    (ixSEG%H06_frac_Sdead  )%varName =trim(cData)   ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
   case('<varname_H06_E_rel_ini>'  ); meta_SEG    (ixSEG%H06_E_rel_ini   )%varName =trim(cData)   ! Hanasaki 2006; Initial release coefficient [-]
   case('<varname_H06_I_Jan>'      ); meta_SEG    (ixSEG%H06_I_Jan       )%varName =trim(cData)   ! Hanasaki 2006; Average January   inflow [m3/s]
   case('<varname_H06_I_Feb>'      ); meta_SEG    (ixSEG%H06_I_Feb       )%varName =trim(cData)   ! Hanasaki 2006; Average Februrary inflow [m3/s]
   case('<varname_H06_I_Mar>'      ); meta_SEG    (ixSEG%H06_I_Mar       )%varName =trim(cData)   ! Hanasaki 2006; Average March     inflow [m3/s]
   case('<varname_H06_I_Apr>'      ); meta_SEG    (ixSEG%H06_I_Apr       )%varName =trim(cData)   ! Hanasaki 2006; Average April     inflow [m3/s]
   case('<varname_H06_I_May>'      ); meta_SEG    (ixSEG%H06_I_May       )%varName =trim(cData)   ! Hanasaki 2006; Average May       inflow [m3/s]
   case('<varname_H06_I_Jun>'      ); meta_SEG    (ixSEG%H06_I_Jun       )%varName =trim(cData)   ! Hanasaki 2006; Average June      inflow [m3/s]
   case('<varname_H06_I_Jul>'      ); meta_SEG    (ixSEG%H06_I_Jul       )%varName =trim(cData)   ! Hanasaki 2006; Average July      inflow [m3/s]
   case('<varname_H06_I_Aug>'      ); meta_SEG    (ixSEG%H06_I_Aug       )%varName =trim(cData)   ! Hanasaki 2006; Average August    inflow [m3/s]
   case('<varname_H06_I_Sep>'      ); meta_SEG    (ixSEG%H06_I_Sep       )%varName =trim(cData)   ! Hanasaki 2006; Average September inflow [m3/s]
   case('<varname_H06_I_Oct>'      ); meta_SEG    (ixSEG%H06_I_Oct       )%varName =trim(cData)   ! Hanasaki 2006; Average October   inflow [m3/s]
   case('<varname_H06_I_Nov>'      ); meta_SEG    (ixSEG%H06_I_Nov       )%varName =trim(cData)   ! Hanasaki 2006; Average November  inflow [m3/s]
   case('<varname_H06_I_Dec>'      ); meta_SEG    (ixSEG%H06_I_Dec       )%varName =trim(cData)   ! Hanasaki 2006; Average December  inflow [m3/s]
   case('<varname_H06_D_Jan>'      ); meta_SEG    (ixSEG%H06_D_Jan       )%varName =trim(cData)   ! Hanasaki 2006; Average January   demand [m3/s]
   case('<varname_H06_D_Feb>'      ); meta_SEG    (ixSEG%H06_D_Feb       )%varName =trim(cData)   ! Hanasaki 2006; Average Februrary demand [m3/s]
   case('<varname_H06_D_Mar>'      ); meta_SEG    (ixSEG%H06_D_Mar       )%varName =trim(cData)   ! Hanasaki 2006; Average March     demand [m3/s]
   case('<varname_H06_D_Apr>'      ); meta_SEG    (ixSEG%H06_D_Apr       )%varName =trim(cData)   ! Hanasaki 2006; Average April     demand [m3/s]
   case('<varname_H06_D_May>'      ); meta_SEG    (ixSEG%H06_D_May       )%varName =trim(cData)   ! Hanasaki 2006; Average May       demand [m3/s]
   case('<varname_H06_D_Jun>'      ); meta_SEG    (ixSEG%H06_D_Jun       )%varName =trim(cData)   ! Hanasaki 2006; Average June      demand [m3/s]
   case('<varname_H06_D_Jul>'      ); meta_SEG    (ixSEG%H06_D_Jul       )%varName =trim(cData)   ! Hanasaki 2006; Average July      demand [m3/s]
   case('<varname_H06_D_Aug>'      ); meta_SEG    (ixSEG%H06_D_Aug       )%varName =trim(cData)   ! Hanasaki 2006; Average August    demand [m3/s]
   case('<varname_H06_D_Sep>'      ); meta_SEG    (ixSEG%H06_D_Sep       )%varName =trim(cData)   ! Hanasaki 2006; Average September demand [m3/s]
   case('<varname_H06_D_Oct>'      ); meta_SEG    (ixSEG%H06_D_Oct       )%varName =trim(cData)   ! Hanasaki 2006; Average October   demand [m3/s]
   case('<varname_H06_D_Nov>'      ); meta_SEG    (ixSEG%H06_D_Nov       )%varName =trim(cData)   ! Hanasaki 2006; Average November  demand [m3/s]
   case('<varname_H06_D_Dec>'      ); meta_SEG    (ixSEG%H06_D_Dec       )%varName =trim(cData)   ! Hanasaki 2006; Average December  demand [m3/s]
   case('<varname_H06_purpose>'    ); meta_SEG    (ixSEG%H06_purpose     )%varName =trim(cData)   ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
   case('<varname_H06_I_mem_F>'    ); meta_SEG    (ixSEG%H06_I_mem_F     )%varName =trim(cData)   ! Hanasaki 2006; Flag to transition to modelled inflow [-]
   case('<varname_H06_D_mem_F>'    ); meta_SEG    (ixSEG%H06_D_mem_F     )%varName =trim(cData)   ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
   case('<varname_H06_I_mem_L>'    ); meta_SEG    (ixSEG%H06_I_mem_L     )%varName =trim(cData)   ! Hanasaki 2006; Memory length in years for inflow [year]
   case('<varname_H06_D_mem_L>'    ); meta_SEG    (ixSEG%H06_D_mem_L     )%varName =trim(cData)   ! Hanasaki 2006; Memory length in years for demand [year]

   ! network topology
   case('<varname_hruContribIx>' ); meta_NTOPO  (ixNTOPO%hruContribIx  )%varName =trim(cData)  ! indices of the vector of HRUs that contribute flow to each segment
   case('<varname_hruContribId>' ); meta_NTOPO  (ixNTOPO%hruContribId  )%varName =trim(cData)  ! ids of the vector of HRUs that contribute flow to each segment
   case('<varname_segId>'        ); meta_NTOPO  (ixNTOPO%segId         )%varName =trim(cData)  ! unique id of each stream segment
   case('<varname_segIndex>'     ); meta_NTOPO  (ixNTOPO%segIndex      )%varName =trim(cData)  ! index of each stream segment (1, 2, 3, ..., n)
   case('<varname_downSegId>'    ); meta_NTOPO  (ixNTOPO%downSegId     )%varName =trim(cData)  ! unique id of the next downstream segment
   case('<varname_downSegIndex>' ); meta_NTOPO  (ixNTOPO%downSegIndex  )%varName =trim(cData)  ! index of downstream reach index
   case('<varname_upSegIds>'     ); meta_NTOPO  (ixNTOPO%upSegIds      )%varName =trim(cData)  ! ids for the vector of immediate upstream stream segments
   case('<varname_upSegIndices>' ); meta_NTOPO  (ixNTOPO%upSegIndices  )%varName =trim(cData)  ! indices for the vector of immediate upstream stream segments
   case('<varname_rchOrder>'     ); meta_NTOPO  (ixNTOPO%rchOrder      )%varName =trim(cData)  ! order that stream segments are processed`
   case('<varname_lakeId>'       ); meta_NTOPO  (ixNTOPO%lakeId        )%varName =trim(cData)  ! unique id of each lake in the river network
   case('<varname_lakeIndex>'    ); meta_NTOPO  (ixNTOPO%lakeIndex     )%varName =trim(cData)  ! index of each lake in the river network
   case('<varname_isLakeInlet>'  ); meta_NTOPO  (ixNTOPO%isLakeInlet   )%varName =trim(cData)  ! flag to define if reach is a lake inlet (1=inlet, 0 otherwise)
   case('<varname_islake>'       ); meta_NTOPO  (ixNTOPO%islake        )%varName =trim(cData)  ! flag to define a lake (1=lake, 0=reach)
   case('<varname_lakeModelType>'); meta_NTOPO  (ixNTOPO%lakeModelType )%varName =trim(cData)  ! defines the lake model type (1=Doll, 2=Hanasaki, 3=etc)
   case('<varname_LakeTargVol>'  ); meta_NTOPO  (ixNTOPO%LakeTargVol   )%varName =trim(cData)  ! flag to follow the provided target volume (1=yes, 0=no)
   case('<varname_userTake>'     ); meta_NTOPO  (ixNTOPO%userTake      )%varName =trim(cData)  ! flag to define if user takes water from reach (1=extract, 0 otherwise)
   case('<varname_goodBasin>'    ); meta_NTOPO  (ixNTOPO%goodBasin     )%varName =trim(cData)  ! flag to define a good basin (1=good, 0=bad)

   ! pfafstetter code
   case('<varname_pfafCode>'     ); meta_PFAF   (ixPFAF%code           )%varName =trim(cData)  ! pfafstetter code

   ! if not in list then keep going
   case default
    message=trim(message)//'unexpected text in control file provided: '//trim(cName)&
                         //' (note strings in control file must match the variable names in public_var.f90)'
    err=20; return
  end select

  ! check I/O error
  if(io_error/=0)then
    message=trim(message)//'problem with internal read of '//trim(cName); err=20; return
  endif

 end do  ! looping through lines in the control file

 ! ---------- directory option  ---------------------------------------------------------------------
 if (trim(restart_dir)==charMissing) then
   restart_dir = output_dir
 endif

 ! ---------- control river network writing option  ---------------------------------------------------------------------
 ! Case1- river network subset mode (idSegOut>0):  Write the network variables read from file over only upstream network specified idSegOut
 ! Case2- river network augment mode: Write full network variables over the entire network
 ! River network subset mode turnes off augmentation mode.
 ! Turned off ntopAugmentMode
 if (idSegOut>0) then
   ntopAugmentMode = .false.
 endif

 ! ---------- time variables  --------------------------------------------------------------------------------------------
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- calendar --- '
   if (trim(calendar)/=charMissing) then
     write(iulog,'(a)') '  calendar is provided in control file: '//trim(calendar)
   else
     write(iulog,'(a)') '  calendar will be read from '//trim(fname_qsim)
   end if
   if (trim(time_units)/=charMissing) then
     write(iulog,'(a)') '  time_unit is provided in control file: '//trim(time_units)
   else
     write(iulog,'(a)') '  time_unit will be read from '//trim(fname_qsim)
   end if
 end if

 ! ---------- runoff unit conversion --------------------------------------------------------------------------------------------
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- runoff unit --- '
   write(iulog,'(a)') '  runoff unit is provided as: '//trim(units_qsim)
 end if
 ! find the position of the "/" character
 ipos = index(trim(units_qsim),'/')
 if(ipos==0)then
   message=trim(message)//'expect the character "/" exists in the units string [units='//trim(units_qsim)//']'
   err=80; return
 endif

 ! get the length and time units
 cLength = units_qsim(1:ipos-1)
 cTime   = units_qsim(ipos+1:len_trim(units_qsim))

 ! get the conversion factor for length
 select case(trim(cLength))
   case('m');  length_conv = 1._dp
   case('mm'); length_conv = 1._dp/1000._dp
   case default
     message=trim(message)//'expect the length units of runoff to be "m" or "mm" [units='//trim(cLength)//']'
     err=81; return
 end select

 ! get the conversion factor for time
 select case(trim(cTime))
   case('d','day');          time_conv = 1._dp/secprday
   case('h','hr','hour');    time_conv = 1._dp/secprhour
   case('s','sec','second'); time_conv = 1._dp
   case default
     message=trim(message)//'expect the time units of to be "day"("d"), "hour"("h") or "second"("s") [time units = '//trim(cTime)//']'
     err=81; return
 end select

 ! ---------- output options --------------------------------------------------------------------------------------------
 ! Assign index for each active routing method
 ! Make sure to turn off write option for routines not used
 if (trim(routOpt)=='0')then; write(iulog,'(a)') 'WARNING: routOpt=0 is accumRunoff option now. 12 is previous 0 now'; endif
 call char2int(trim(routOpt), routeMethods, invalid_value=0)
 nRoutes = size(routeMethods)
 onRoute = .false.
 do iRoute = 1, nRoutes
   select case(routeMethods(iRoute))
     case(accumRunoff);           idxSUM = iRoute; onRoute(accumRunoff)=.true.
     case(kinematicWaveTracking); idxKWT = iRoute; onRoute(kinematicWaveTracking)=.true.
     case(impulseResponseFunc);   idxIRF = iRoute; onRoute(impulseResponseFunc)=.true.
     case(muskingumCunge);        idxMC  = iRoute; onRoute(muskingumCunge)=.true.
     case(kinematicWave);         idxKW  = iRoute; onRoute(kinematicWave)=.true.
     case(diffusiveWave);         idxDW  = iRoute; onRoute(diffusiveWave)=.true.
     case default
       message=trim(message)//'routOpt may include invalid digits; expect digits 1-5 in routOpt'; err=81; return
   end select
 end do

 do iRoute = 0, nRouteMethods-1
   select case(iRoute)
     case(accumRunoff);           if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile=.false.
     case(kinematicWaveTracking); if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%KWTroutedRunoff)%varFile=.false.
     case(impulseResponseFunc);   if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%IRFroutedRunoff)%varFile=.false.
     case(muskingumCunge);        if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%MCroutedRunoff)%varFile=.false.
     case(kinematicWave);         if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%KWroutedRunoff)%varFile=.false.
     case(diffusiveWave);         if (.not. onRoute(iRoute)) meta_rflx(ixRFLX%DWroutedRunoff)%varFile=.false.
     case default; message=trim(message)//'expect digits from 1 and 5'; err=81; return
   end select
 end do

 ! basin runoff routing option
 if (doesBasinRoute==0) meta_rflx(ixRFLX%instRunoff)%varFile = .false.

 END SUBROUTINE read_control

END MODULE read_control_module
