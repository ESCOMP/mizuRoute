module read_control_module
USE nrtype
USE public_var
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
 USE public_var
 USE globalData, only:time_conv,length_conv  ! conversion factors

 ! metadata structures
 USE globalData, only : meta_HRU             ! HRU properties
 USE globalData, only : meta_HRU2SEG         ! HRU-to-segment mapping
 USE globalData, only : meta_SEG             ! stream segment properties
 USE globalData, only : meta_NTOPO           ! network topology

 ! named variables in each structure
 USE var_lookup, only : ixHRU                ! index of variables for data structure
 USE var_lookup, only : ixHRU2SEG            ! index of variables for data structure
 USE var_lookup, only : ixSEG                ! index of variables for data structure
 USE var_lookup, only : ixNTOPO              ! index of variables for data structure

 ! external subroutines
 USE ascii_util_module,only:file_open        ! open file (performs a few checks as well)
 USE ascii_util_module,only:get_vlines       ! get a list of character strings from non-comment lines

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
   print*, trim(cName), ' --> ', trim(cData)

   ! populate variables
   select case(trim(cName))

   ! part 1: define directories
   case('<ancil_dir>');         ancil_dir                          = trim(cData) ! directory containing ancillary data
   case('<input_dir>');         input_dir                          = trim(cData) ! directory containing input data
   case('<output_dir>');        output_dir                         = trim(cData) ! directory containing output data

   ! part 2: Define river network netCDF and its meta
   case('<fname_ntop>');        fname_ntop                               = trim(cData) ! name of file containing stream network topology information
   case('<fname_sseg>');        fname_sseg                               = trim(cData) ! name of file containing stream segment information
   case('<dname_nhru>');        dname_nhru                               = trim(cData) ! dimension name of the HRUs
   case('<dname_sseg>');        dname_sseg                               = trim(cData) ! dimension name of the stream segments
   ! part 3: Define desired outlet segment ID (if -9999 --> route over the entire network)
   case('<seg_outlet>'   );     read(cData,*,iostat=io_error) idSegOut           ! seg_id of outlet streamflow segment

   ! Part 4: Define runoff netCDF and its meta
   case('<fname_qsim>');        fname_qsim   = trim(cData)                       ! name of file containing the runoff
   case('<vname_qsim>');        vname_qsim   = trim(cData)                       ! name of runoff variable
   case('<vname_time>');        vname_time   = trim(cData)                       ! name of time variable in the runoff file
   case('<vname_hruid>');       vname_hruid  = trim(cData)                       ! name of the HRU id
   case('<dname_time>');        dname_time   = trim(cData)                       ! name of time variable in the runoff file
   case('<dname_hruid>');       dname_hruid  = trim(cData)                       ! name of the HRU id
   case('<units_qsim>');        units_qsim   = trim(cData)                       ! units of runoff
   case('<dt_qsim>');           read(cData,*,iostat=io_error) dt                 ! time interval of the gridded runoff

   ! Part 5: Define runoff mapping netCDF file and its metadata
   case('<is_remap>');             read(cData,*,iostat=io_error) is_remap        ! logical whether or not runnoff needs to be mapped to river network HRU
   case('<fname_remap>');          fname_remap           = trim(cData)           ! name of runoff mapping netCDF
   case('<vname_hruid_in_remap>'); vname_hruid_in_remap  = trim(cData)           ! name of variable containing ID of runoff HRU
   case('<vname_weight>');         vname_weight          = trim(cData)           ! name of variable contating areal weights of runoff HRUs within each river network
   case('<vname_qhruid>');         vname_qhruid          = trim(cData)           ! name of variable containing ID of runoff HRU
   case('<vname_num_qhru>');       vname_num_qhru        = trim(cData)           ! name of variable containing numbers of runoff HRUs within each river network HRU
   case('<dname_hru_remap>');      dname_hru_remap       = trim(cData)           ! name of variable containing ID of runoff HRU
   case('<dname_data_remap>');     dname_data_remap      = trim(cData)           ! name of variable containing ID of runoff HRU

   ! part 6 Define run control
   case('<restart_opt>');   read(cData,*,iostat=io_error) isRestart              ! restart option: True-> model run with restart, F -> model run with empty channels
   case('<route_opt>');     read(cData,*,iostat=io_error) routOpt                ! routing scheme options  0-> both, 1->IRF, 2->KWT, otherwise error

   ! Part 7: Define output filename
   case('<fname_output>');    fname_output    = trim(cData)                      ! filename for the model output
   case('<fname_state_in>');  fname_state_in  = trim(cData)                      ! filename for the channel states
   case('<fname_state_out>'); fname_state_out = trim(cData)                      ! filename for the channel states

   ! Part 8: Define options to include/skip calculations
   case('<hydGeometryOption>'); read(cData,*,iostat=io_error) hydGeometryOption  ! option for hydraulic geometry calculations (0=read from file, 1=compute)
   case('<topoNetworkOption>'); read(cData,*,iostat=io_error) topoNetworkOption  ! option for network topology calculations (0=read from file, 1=compute)
   case('<computeReachList>');  read(cData,*,iostat=io_error) computeReachList   ! option to compute list of upstream reaches (0=do not compute, 1=compute)

   ! Part 9: Define namelist name for routing parameters
   case('<param_nml>');       param_nml       = trim(cData)                      ! name of namelist including routing parameter value

   ! Part X: Dictionary to modify variable names

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
   case('<varname_upsLength>'    ); meta_SEG    (ixSEG%upsLength       )%varName =trim(cData)  ! length of the vector of reaches above each reach (m)
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

   ! if not in list then keep going
   case default
    message=trim(message)//'unexpected text in control file -- provided '//trim(cName)&
                         //' (note strings in control file must match the variable names in var_lookup.f90)'
    err=20; return

  end select

  ! check I/O error
  if(io_error/=0)then
   message=trim(message)//'problem with internal read of '//trim(cName)
   err=20; return
  endif

 end do  ! looping through lines in the control file

 ! ---------- unit conversion --------------------------------------------------------------------------------------------

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
  case('d','day');    time_conv = 1._dp/secprday
  case('h','hour');   time_conv = 1._dp/secprhour
  case('s','second'); time_conv = 1._dp
  case default
   message=trim(message)//'cannot identify the time units [time units = '//trim(cTime)//']'
   err=81; return
 end select

 end subroutine read_control

end module read_control_module
