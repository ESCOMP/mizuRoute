module read_control_module

USE nrtype
USE public_var
USE globalData, ONLY: pid, nNodes, masterproc             ! procs id and number of procs

implicit none
! privacy
private
public::read_control
contains

 ! =======================================================================================================
 ! * new subroutine: read the control file
 ! =======================================================================================================
 ! read the control file
 subroutine read_control(ctl_fname, err, message)

 ! data types
 USE nrtype                                  ! variable types, etc.

 ! global vars
 USE globalData, ONLY: time_conv,length_conv  ! conversion factors

 ! metadata structures
 USE globalData, ONLY: meta_HRU             ! HRU properties
 USE globalData, ONLY: meta_HRU2SEG         ! HRU-to-segment mapping
 USE globalData, ONLY: meta_SEG             ! stream segment properties
 USE globalData, ONLY: meta_NTOPO           ! network topology
 USE globalData, ONLY: meta_PFAF            ! pfafstetter code
 USE globalData, ONLY: meta_rflx            ! river flux variables

 ! named variables in each structure
 USE var_lookup, ONLY: ixHRU                ! index of variables for data structure
 USE var_lookup, ONLY: ixHRU2SEG            ! index of variables for data structure
 USE var_lookup, ONLY: ixSEG                ! index of variables for data structure
 USE var_lookup, ONLY: ixNTOPO              ! index of variables for data structure
 USE var_lookup, ONLY: ixPFAF               ! index of variables for data structure
 USE var_lookup, ONLY: ixRFLX               ! index of variables for data structure

 ! external subroutines
 USE ascii_util_module, ONLY:file_open        ! open file (performs a few checks as well)
 USE ascii_util_module, ONLY:get_vlines       ! get a list of character strings from non-comment lines

 implicit none
 ! input
 character(*), intent(in)          :: ctl_fname      ! name of the control file
 ! output: error control
 integer(i4b),intent(out)          :: err            ! error code
 character(*),intent(out)          :: message        ! error message
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
 character(len=strLen)             :: cmessage       ! error message from subroutine
 ! initialize error control
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
 if (masterproc) then
   write(iulog,'(2a)') new_line('a'), '---- read control file --- '
 end if
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
   case('<ancil_dir>');            ancil_dir   = trim(cData)                       ! directory containing ancillary data
   case('<input_dir>');            input_dir   = trim(cData)                       ! directory containing input data
   case('<output_dir>');           output_dir  = trim(cData)                       ! directory containing output data
   ! RUN CONTROL
   case('<case_name>');            case_name   = trim(cData)                       ! name of simulation. used as head of model output and restart file
   case('<sim_start>');            simStart    = trim(cData)                       ! date string defining the start of the simulation
   case('<sim_end>');              simEnd      = trim(cData)                       ! date string defining the end of the simulation
   case('<newFileFrequency>');     newFileFrequency = trim(cData)                  ! frequency for new output files (day, month, annual, single)
   case('<restart_write>');        restart_write        = trim(cData)              ! restart write option: N[n]ever, L[l]ast
   case('<restart_date>');         restart_date         = trim(cData)              ! specified restart date, yyyy-mm-dd (hh:mm:ss)
   case('<fname_state_in>');       fname_state_in  = trim(cData)                   ! filename for the channel states
   case('<route_opt>');            read(cData,*,iostat=io_error) routOpt           ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error
   case('<doesBasinRoute>');       read(cData,*,iostat=io_error) doesBasinRoute    ! basin routing options   0-> no, 1->IRF, otherwise error
   case('<doesAccumRunoff>');      read(cData,*,iostat=io_error) doesAccumRunoff   ! option to delayed runoff accumulation over all the upstream reaches. 0->no, 1->yes
   case('<seg_outlet>'   );        read(cData,*,iostat=io_error) idSegOut          ! desired outlet reach id (if -9999 --> route over the entire network)
   ! RIVER NETWORK TOPOLOGY
   case('<fname_ntopOld>');        fname_ntopOld = trim(cData)                     ! name of file containing stream network topology information
   case('<ntopAugmentMode>');      read(cData,*,iostat=io_error) ntopAugmentMode   ! option for river network augmentation mode. terminate the program after writing augmented ntopo.
   case('<fname_ntopNew>');        fname_ntopNew = trim(cData)                     ! name of file containing stream segment information
   case('<dname_nhru>');           dname_nhru    = trim(cData)                     ! dimension name of the HRUs
   case('<dname_sseg>');           dname_sseg    = trim(cData)                     ! dimension name of the stream segments
   ! RUNOFF FILE
   case('<fname_qsim>');           fname_qsim   = trim(cData)                      ! name of text file containing nc file names and their order for runoff
   case('<vname_qsim>');           vname_qsim   = trim(cData)                      ! name of runoff variable
   case('<vname_time>');           vname_time   = trim(cData)                      ! name of time variable in the runoff file
   case('<vname_hruid>');          vname_hruid  = trim(cData)                      ! name of the HRU id
   case('<dname_time>');           dname_time   = trim(cData)                      ! name of time variable in the runoff file
   case('<dname_hruid>');          dname_hruid  = trim(cData)                      ! name of the HRU id dimension
   case('<dname_xlon>');           dname_xlon   = trim(cData)                      ! name of x (j,lon) dimension
   case('<dname_ylat>');           dname_ylat   = trim(cData)                      ! name of y (i,lat) dimension
   case('<units_qsim>');           units_qsim   = trim(cData)                      ! units of runoff
   case('<dt_qsim>');              read(cData,*,iostat=io_error) dt                ! time interval of the gridded runoff
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
   ! SPATIAL CONSTANT PARAMETERS
   case('<param_nml>');            param_nml       = trim(cData)                   ! name of namelist including routing parameter value
   ! USER OPTIONS: Define options to include/skip calculations
   case('<hydGeometryOption>');    read(cData,*,iostat=io_error) hydGeometryOption ! option for hydraulic geometry calculations (0=read from file, 1=compute)
   case('<topoNetworkOption>');    read(cData,*,iostat=io_error) topoNetworkOption ! option for network topology calculations (0=read from file, 1=compute)
   case('<computeReachList>');     read(cData,*,iostat=io_error) computeReachList  ! option to compute list of upstream reaches (0=do not compute, 1=compute)
   ! TIME
   case('<time_units>');           time_units = trim(cData)                        ! time units. format should be <unit> since yyyy-mm-dd (hh:mm:ss). () can be omitted
   case('<calendar>');             calendar   = trim(cData)                        ! calendar name
   ! MISCELLANEOUS
   case('<debug>');                read(cData,*,iostat=io_error) debug             ! print out detailed information throught the probram
   case('<desireId>'   );          read(cData,*,iostat=io_error) desireId          ! turn off checks or speficy reach ID if necessary to print on screen
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
   message=trim(message)//'problem with internal read of '//trim(cName)
   err=20; return
  endif

 end do  ! looping through lines in the control file

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

 end subroutine read_control

end module read_control_module
