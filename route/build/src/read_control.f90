MODULE read_control_module

! <overall comments>

USE nrtype
USE public_var

implicit none

private
public::read_control

CONTAINS

 ! =======================================================================================================
 ! subroutine: read the control file
 ! =======================================================================================================
 SUBROUTINE read_control(ctl_fname, err, message)

 ! global vars
 USE globalData, only:time_conv,length_conv  ! conversion factors

 ! metadata structures
 USE globalData, ONLY: meta_HRU             ! HRU properties
 USE globalData, ONLY: meta_HRU2SEG         ! HRU-to-segment mapping
 USE globalData, ONLY: meta_SEG             ! stream segment properties
 USE globalData, ONLY: meta_NTOPO           ! network topology
 USE globalData, ONLY: meta_PFAF            ! pfafstetter code
 USE globalData, ONLY: meta_rflx            ! river flux variables
 USE globalData, ONLY: nRoutes              ! number of active routing methods
 USE globalData, ONLY: routeMethods         ! active routing method index and id
 USE globalData, ONLY: onRoute                 ! logical to indicate actiive routing method(s)
 USE globalData, ONLY: idxSUM,idxIRF,idxKWT, &
                        idxKW, idxMC, idxDW
 ! named variables in each structure
 USE var_lookup, ONLY : ixHRU                ! index of variables for data structure
 USE var_lookup, ONLY : ixHRU2SEG            ! index of variables for data structure
 USE var_lookup, ONLY : ixSEG                ! index of variables for data structure
 USE var_lookup, ONLY : ixNTOPO              ! index of variables for data structure
 USE var_lookup, ONLY : ixPFAF               ! index of variables for data structure
 USE var_lookup, ONLY : ixRFLX, nVarsRFLX    ! index of variables for data structure

 ! external subroutines
 USE ascii_util_module, only: file_open      ! open file (performs a few checks as well)
 USE ascii_util_module, only: get_vlines     ! get a list of character strings from non-comment lines
 USE nr_utility_module, ONLY: char2int       ! convert integer number to a array containing individual digits

 implicit none
 ! arguments
 character(*), intent(in)          :: ctl_fname      ! name of the control file
 integer(i4b), intent(out)         :: err            ! error code
 character(*), intent(out)         :: message        ! error message
 ! ------------------------------------------------------------------------------------------------------
 ! Local variables
 character(len=strLen),allocatable :: cLines(:)      ! vector of character strings
 character(len=strLen)             :: cName,cData    ! name and data from cLines(iLine)
 character(len=strLen)             :: cLength,cTime  ! length and time units
 integer(i4b)                      :: ipos           ! index of character string
 integer(i4b)                      :: ibeg_name      ! start index of variable name in string cLines(iLine)
 integer(i4b)                      :: iend_name      ! end index of variable name in string cLines(iLine)
 integer(i4b)                      :: iend_data      ! end index of data in string cLines(iLine)
 integer(i4b)                      :: iLine          ! index of line in cLines
 integer(i4b)                      :: iunit          ! file unit
 integer(i4b)                      :: io_error       ! error in I/O
 integer(i4b)                      :: iRoute         ! loop index
 character(len=strLen)             :: cmessage       ! error message from subroutine

 err=0; message='read_control/'

 ! *** get a list of character strings from non-comment lines ****
 ! open file (also returns un-used file unit used to open the file)
 call file_open(trim(ctl_fname),iunit,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

 ! get a list of character strings from non-comment lines
 call get_vlines(iunit,cLines,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage);return;endif

 ! close the file unit
 close(iunit)

 ! loop through the non-comment lines in the input file, and extract the name and the information
 write(iulog,'(a)') '---- read control file --- '
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
   write(iulog,'(x,a,a,a)') trim(cName), ' --> ', trim(cData)

   ! populate variables
   select case(trim(cName))

   ! DIRECTORIES
   case('<ancil_dir>');            ancil_dir   = trim(cData)                       ! directory containing ancillary data (network, mapping, namelist)
   case('<input_dir>');            input_dir   = trim(cData)                       ! directory containing input runoff netCDF
   case('<output_dir>');           output_dir  = trim(cData)                       ! directory for routed flow output (netCDF)
   case('<restart_dir>');          restart_dir = trim(cData)                       ! directory for restart output (netCDF)
   ! RIVER NETWORK TOPOLOGY
   case('<fname_ntopOld>');        fname_ntopOld = trim(cData)                     ! name of file containing stream network topology information
   case('<ntopAugmentMode>');      read(cData,*,iostat=io_error) ntopAugmentMode   ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
   case('<fname_ntopNew>');        fname_ntopNew = trim(cData)                     ! name of file containing stream segment information
   case('<dname_nhru>');           dname_nhru    = trim(cData)                     ! dimension name of the HRUs
   case('<dname_sseg>');           dname_sseg    = trim(cData)                     ! dimension name of the stream segments
   ! RUNOFF FILE
   case('<fname_qsim>');           fname_qsim   = trim(cData)                      ! name of file containing the runoff
   case('<vname_qsim>');           vname_qsim   = trim(cData)                      ! name of runoff variable
   case('<vname_time>');           vname_time   = trim(cData)                      ! name of time variable in the runoff file
   case('<vname_hruid>');          vname_hruid  = trim(cData)                      ! name of the HRU id
   case('<dname_time>');           dname_time   = trim(cData)                      ! name of time variable in the runoff file
   case('<dname_hruid>');          dname_hruid  = trim(cData)                      ! name of the HRU id dimension
   case('<dname_xlon>');           dname_xlon   = trim(cData)                      ! name of x (j,lon) dimension
   case('<dname_ylat>');           dname_ylat   = trim(cData)                      ! name of y (i,lat) dimension
   case('<units_qsim>');           units_qsim   = trim(cData)                      ! units of runoff
   case('<dt_qsim>');              read(cData,*,iostat=io_error) dt                ! time interval of the gridded runoff
   case('<ro_fillvalue>')
                                   read(cData,*,iostat=io_error) ro_fillvalue      ! fillvalue used for runoff depth variable
                                   userRunoffFillvalue = .true.                    ! true -> runoff depth fillvalue used in netcdf is specified here, otherwise -> false
   ! RUNOFF REMAPPING
   case('<is_remap>');             read(cData,*,iostat=io_error) is_remap          ! logical whether or not runnoff needs to be mapped to river network HRU
   case('<fname_remap>');          fname_remap          = trim(cData)              ! name of runoff mapping netCDF
   case('<vname_hruid_in_remap>'); vname_hruid_in_remap = trim(cData)              ! name of variable containing ID of river network HRU
   case('<vname_weight>');         vname_weight         = trim(cData)              ! name of variable contating areal weights of runoff HRUs within each river network HRU
   case('<vname_qhruid>');         vname_qhruid         = trim(cData)              ! name of variable containing ID of runoff HRU
   case('<vname_num_qhru>');       vname_num_qhru       = trim(cData)              ! name of variable containing numbers of runoff HRUs within each river network HRU
   case('<vname_i_index>');        vname_i_index        = trim(cData)              ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
   case('<vname_j_index>');        vname_j_index        = trim(cData)              ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)
   case('<dname_hru_remap>');      dname_hru_remap      = trim(cData)              ! name of dimension of river network HRU ID
   case('<dname_data_remap>');     dname_data_remap     = trim(cData)              ! name of dimension of runoff HRU overlapping with river network HRU
   ! RUN CONTROL
   case('<case_name>');            case_name            = trim(cData)              ! name of simulation. used as head of model output and restart file
   case('<sim_start>');            simStart    = trim(cData)                       ! date string defining the start of the simulation
   case('<sim_end>');              simEnd      = trim(cData)                       ! date string defining the end of the simulation
   case('<route_opt>');            routOpt     = trim(cData)                       ! routing scheme options  0-> accumRunoff, 1->IRF, 2->KWT, 3-> KW, 4->MC, 5->DW
   case('<doesBasinRoute>');       read(cData,*,iostat=io_error) doesBasinRoute    ! basin routing options   0-> no, 1->IRF, otherwise error
   case('<newFileFrequency>');     newFileFrequency     = trim(cData)              ! frequency for new output options (case-insensitive): daily, monthly, yearly, or single
   ! RESTART
   case('<restart_write>');        restart_write        = trim(cData)              ! restart write option (case-insensitive): never, last, specified, yearly, monthly, or daily
   case('<restart_date>');         restart_date         = trim(cData)              ! specified restart date, yyyy-mm-dd (hh:mm:ss) for Specified option
   case('<restart_month>');        read(cData,*,iostat=io_error) restart_month     ! restart periodic month
   case('<restart_day>');          read(cData,*,iostat=io_error) restart_day       ! restart periodic day
   case('<restart_hour>');         read(cData,*,iostat=io_error) restart_hour      ! restart periodic hour
   case('<fname_state_in>');       fname_state_in       = trim(cData)              ! filename for the channel states
   ! SPATIAL CONSTANT PARAMETERS
   case('<param_nml>');            param_nml       = trim(cData)                   ! name of namelist including routing parameter value
   ! USER OPTIONS: Define options to include/skip calculations
   case('<qtakeOption>');          read(cData,*,iostat=io_error) qtakeOption       ! option for abstraction/injection option
   case('<hydGeometryOption>');    read(cData,*,iostat=io_error) hydGeometryOption ! option for hydraulic geometry calculations (0=read from file, 1=compute)
   case('<topoNetworkOption>');    read(cData,*,iostat=io_error) topoNetworkOption ! option for network topology calculations (0=read from file, 1=compute)
   case('<computeReachList>');     read(cData,*,iostat=io_error) computeReachList  ! option to compute list of upstream reaches (0=do not compute, 1=compute)
   ! TIME
   case('<time_units>');           time_units = trim(cData)                        ! time units. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
   case('<calendar>');             calendar   = trim(cData)                        ! calendar name
   ! MISCELLANEOUS
   case('<debug>');                read(cData,*,iostat=io_error) debug             ! print out detailed information throught the probram
   case('<seg_outlet>'   );        read(cData,*,iostat=io_error) idSegOut          ! desired outlet reach id (if -9999 --> route over the entire network)
   case('<desireId>'   );          read(cData,*,iostat=io_error) desireId          ! turn off checks or speficy reach ID if necessary to print on screen
   case('<netcdf_format>');        netcdf_format = trim(cData)                     ! netcdf format for output 'classic','64bit_offset','netcdf4'
   ! PFAFCODE
   case('<maxPfafLen>');           read(cData,*,iostat=io_error) maxPfafLen        ! maximum digit of pfafstetter code (default 32)
   case('<pfafMissing>');          pfafMissing = trim(cData)                       ! missing pfafcode (e.g., reach without any upstream area)
   ! OUTPUT OPTIONS
   case('<basRunoff>');            read(cData,*,iostat=io_error) meta_rflx(ixRFLX%basRunoff        )%varFile
   case('<instRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%instRunoff       )%varFile
   case('<dlayRunoff>');           read(cData,*,iostat=io_error) meta_rflx(ixRFLX%dlayRunoff       )%varFile
   case('<sumUpstreamRunoff>');    read(cData,*,iostat=io_error) meta_rflx(ixRFLX%sumUpstreamRunoff)%varFile
   case('<KWTroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWTroutedRunoff  )%varFile
   case('<IRFroutedRunoff>');      read(cData,*,iostat=io_error) meta_rflx(ixRFLX%IRFroutedRunoff  )%varFile
   case('<KWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%KWroutedRunoff   )%varFile
   case('<DWroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%DWroutedRunoff   )%varFile
   case('<MCroutedRunoff>');       read(cData,*,iostat=io_error) meta_rflx(ixRFLX%MCroutedRunoff   )%varFile

   ! VARIABLE NAMES for data (overwrite default name in popMeta.f90)
   ! HRU structure
   case('<varname_area>'         ); meta_HRU    (ixHRU%area            )%varName = trim(cData) ! HRU area
   ! Mapping from HRUs to stream segments
   case('<varname_HRUid>'        ); meta_HRU2SEG(ixHRU2SEG%HRUid       )%varName = trim(cData) ! HRU id
   case('<varname_HRUindex>'     ); meta_HRU2SEG(ixHRU2SEG%HRUindex    )%varName = trim(cData) ! HRU index
   case('<varname_hruSegId>'     ); meta_HRU2SEG(ixHRU2SEG%hruSegId    )%varName = trim(cData) ! the stream segment id below each HRU
   case('<varname_hruSegIndex>'  ); meta_HRU2SEG(ixHRU2SEG%hruSegIndex )%varName = trim(cData) ! the stream segment index below each HRU
   ! reach properties
   case('<varname_length>'       ); meta_SEG    (ixSEG%length          )%varName =trim(cData)  ! length of segment  (m)
   case('<varname_slope>'        ); meta_SEG    (ixSEG%slope           )%varName =trim(cData)  ! slope of segment   (-)
   case('<varname_width>'        ); meta_SEG    (ixSEG%width           )%varName =trim(cData)  ! width of segment   (m)
   case('<varname_man_n>'        ); meta_SEG    (ixSEG%man_n           )%varName =trim(cData)  ! Manning's n        (weird units)
   case('<varname_hruArea>'      ); meta_SEG    (ixSEG%hruArea         )%varName =trim(cData)  ! local basin area (m2)
   case('<varname_weight>'       ); meta_SEG    (ixSEG%weight          )%varName =trim(cData)  ! HRU weight
   case('<varname_timeDelayHist>'); meta_SEG    (ixSEG%timeDelayHist   )%varName =trim(cData)  ! time delay histogram for each reach (s)
   case('<varname_upsArea>'      ); meta_SEG    (ixSEG%upsArea         )%varName =trim(cData)  ! area above the top of the reach -- zero if headwater (m2)
   case('<varname_basUnderLake>' ); meta_SEG    (ixSEG%basUnderLake    )%varName =trim(cData)  ! Area of basin under lake  (m2)
   case('<varname_rchUnderLake>' ); meta_SEG    (ixSEG%rchUnderLake    )%varName =trim(cData)  ! Length of reach under lake (m)
   case('<varname_qtake>'        ); meta_SEG    (ixSEG%Qtake           )%varName =trim(cData)  ! abstraction(-)/injection(+) (m3 s-1)
   case('<varname_minFlow>'      ); meta_SEG    (ixSEG%minFlow         )%varName =trim(cData)  ! minimum environmental flow
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
   case('<varname_userTake>'     ); meta_NTOPO  (ixNTOPO%userTake      )%varName =trim(cData)  ! flag to define if user takes water from reach (1=extract, 0 otherwise)
   case('<varname_goodBasin>'    ); meta_NTOPO  (ixNTOPO%goodBasin     )%varName =trim(cData)  ! flag to define a good basin (1=good, 0=bad)
   ! pfafstetter code
   case('<varname_pfafCode>'     ); meta_PFAF   (ixPFAF%code           )%varName =trim(cData)  ! pfafstetter code

   ! if not in list then keep going
   case default
    message=trim(message)//'unexpected text in control file -- provided '//trim(cName)&
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

 ! ---------- runoff unit conversion --------------------------------------------------------------------------------------------

 write(iulog,'(2a)') new_line('a'), '---- runoff unit --- '
 write(iulog,'(a)') '  runoff unit is provided as: '//trim(units_qsim)
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
     message=trim(message)//'expect the length units to be "m" or "mm" [units='//trim(cLength)//']'
     err=81; return
 end select

 ! get the conversion factor for time
 select case(trim(cTime))
   case('d','day');          time_conv = 1._dp/secprday
   case('h','hr','hour');    time_conv = 1._dp/secprhour
   case('s','sec','second'); time_conv = 1._dp
   case default
     message=trim(message)//'expect the time units to be "day"("d"), "hour"("h") or "second"("s") [time units = '//trim(cTime)//']'
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
     case default; message=trim(message)//'expect digits from 0 and 5'; err=81; return
   end select
 end do

 ! basin runoff routing option
 if (doesBasinRoute==0) meta_rflx(ixRFLX%instRunoff)%varFile=.false.

 END SUBROUTINE read_control

END MODULE read_control_module
