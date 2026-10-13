"""
Python module to read in a sample mizuRoute control file into a data type
as an array of keys as well as a dictionary of the values. The Dictionary
can then be modified and output as a new file.

Erik Kluzek
"""

import sys, re, os, logging, collections

sys.path.append( "../../cime/scripts/lib" );
sys.path.append( "../../../../cime/scripts/lib" );

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type, run_cmd_no_fail

try:
    import tomllib
except ImportError:
    import tomli as tomllib

logger = logging.getLogger(__name__)

class mizuRoute_control(object):

   """ Object to hold a dictionary of settings for mizuRoute control """

   # Class Data:
   fileRead = False                         # If file has been read or not
   lineMatch = '^<(.+?)>\s+(.*?)\s*\!(.*)$' # Pattern to match for legacy lines
   longestName = 0                          # Longest name
   longestValue = 0                         # Longest value

   def __init__(self):
      self.ctldict = collections.OrderedDict()     # Ordered dictionary of control elements
      self.comments = {}                           # Comments associated with keys

   @classmethod
   def from_toml( cls, infile, allowEmpty=False ):
       """ Factory constructor to read and parse a TOML format control file """
       inst = cls()
       inst.readToml( infile, allowEmpty=allowEmpty )
       return inst

   @classmethod
   def from_control( cls, infile, allowEmpty=False ):
       """ Factory constructor to read and parse a legacy format control file """
       inst = cls()
       inst.readControl( infile, allowEmpty=allowEmpty )
       return inst

   def read( self, infile, allowEmpty=False ):
       """
       Read and parse a mizuRoute control file (auto-detect format)
       """
       if ( infile.endswith(".toml") or infile.endswith("_toml") or os.path.basename(infile) == "user_nl_mizuroute_toml" ):
           return self.readToml( infile, allowEmpty=allowEmpty )
       else:
           return self.readControl( infile, allowEmpty=allowEmpty )

   def readControl( self, infile, allowEmpty=False ):
       """
       Read and parse a legacy mizuRoute control file
       """
       logger.debug( "read in legacy control file: "+infile )
       if ( not os.path.exists(infile) ):
          expect( False, "Input file to read does NOT exist: "+infile )

       ctlfile = open( infile, "r" )
       lines = ctlfile.readlines()
       ctlfile.close()

       # Loop through each line in the file
       for line in lines:
          # Ignore comment lines
          if ( not line.find( "!" ) == 0 and line.strip() ):
             match = re.search( self.lineMatch, line )
             if ( not match ):
                expect( False, "Error in reading in line:"+line )
             else:
                name = match.group(1).strip()
                value = match.group(2).strip()
                comment = match.group(3).strip() if len(match.groups()) >= 3 and match.group(3) else ""
                self.set( name, value, allowNewName=True, comment=comment )

       # If no data was read -- abort with an error
       if ( len(self.ctldict) == 0 and not allowEmpty ):
          expect( False, "No data was read from the file: "+infile )

       # Mark the file as read
       logger.debug( "File read" )
       self.fileRead = True

   VALID_CONTROL_KEYS = {
        'ancil_dir', 'input_dir', 'output_dir', 'restart_dir', 'case_name',
        'sim_start', 'sim_end', 'continue_run', 'route_opt', 'doesBasinRoute',
        'dt_qsim', 'floodplain', 'hw_drain_point', 'tracer', 'is_lake_sim',
        'lakeRegulate', 'LakeInputOption', 'is_flux_wm', 'is_vol_wm',
        'is_vol_wm_jumpstart', 'scale_factor_runoff', 'offset_value_runoff',
        'scale_factor_Ep', 'offset_value_Ep', 'is_Ep_upward_negative',
        'scale_factor_prec', 'offset_value_prec', 'min_length_route',
        'fname_ntopOld', 'ntopAugmentMode', 'fname_ntopNew', 'dname_nhru', 'dname_sseg',
        'fname_qsim', 'vname_qsim', 'vname_evapo', 'vname_precip', 'vname_solute',
        'vname_time', 'vname_hruid', 'dname_time', 'dname_hruid', 'dname_xlon', 'dname_ylat',
        'units_qsim', 'units_cc', 'dt_ro', 'input_fillvalue', 'ro_calendar',
        'ro_time_units', 'ro_time_stamp', 'runoffMin', 'fname_wm', 'vname_flux_wm',
        'vname_vol_wm', 'vname_time_wm', 'vname_segid_wm', 'dname_time_wm',
        'dname_segid_wm', 'dt_wm', 'is_remap', 'fname_remap', 'vname_hruid_in_remap',
        'vname_weight', 'vname_qhruid', 'vname_num_qhru', 'vname_i_index', 'vname_j_index',
        'dname_hru_remap', 'dname_data_remap', 'restart_write', 'restart_date',
        'restart_month', 'restart_day', 'restart_hour', 'fname_state_in', 'param_nml',
        'qmodOption', 'qBlendPeriod', 'QerrTrend', 'hydGeometryOption', 'topoNetworkOption',
        'computeReachList', 'gageMetaFile', 'outputAtGage', 'fname_gageObs', 'vname_gageFlow',
        'vname_gageSite', 'vname_gageTime', 'dname_gageSite', 'dname_gageTime', 'strlen_gageSite',
        'pio_netcdf_format', 'pio_netcdf_type', 'debug', 'seg_outlet', 'desireId',
        'checkMassBalance', 'maxPfafLen', 'pfafMissing', 'time_units', 'newFileFrequency',
        'outputFrequency', 'outputNameOption', 'histTimeStamp_offset', 'outputInflow',
        'qgwl_runoff_option', 'bypass_routing_option', 'correct_area', 'ice_runoff',
        'varname_area', 'varname_HRUid', 'varname_HRUindex', 'varname_hruSegId',
        'varname_hruSegIndex', 'varname_length', 'varname_slope', 'varname_width',
        'varname_depth', 'varname_sideSlope', 'varname_man_n', 'varname_floodplainSlope',
        'varname_hruArea', 'varname_weight', 'varname_timeDelayHist', 'varname_upsArea',
        'varname_basUnderLake', 'varname_rchUnderLake', 'varname_minFlow', 'varname_D03_MaxStorage',
        'varname_D03_Coefficient', 'varname_D03_Power', 'varname_D03_S0', 'varname_HYP_E_emr',
        'varname_HYP_E_lim', 'varname_HYP_E_min', 'varname_HYP_E_zero', 'varname_HYP_Qrate_emr',
        'varname_HYP_Erate_emr', 'varname_HYP_Qrate_prim', 'varname_HYP_Qrate_amp',
        'varname_HYP_Qrate_phs', 'varname_HYP_prim_F', 'varname_HYP_A_avg', 'varname_HYP_Qsim_mode',
        'varname_H06_Smax', 'varname_H06_alpha', 'varname_H06_envfact', 'varname_H06_S_ini',
        'varname_H06_c1', 'varname_H06_c2', 'varname_H06_exponent', 'varname_H06_denominator',
        'varname_H06_c_compare', 'varname_H06_frac_Sdead', 'varname_H06_E_rel_ini',
        'varname_H06_I_Jan', 'varname_H06_I_Feb', 'varname_H06_I_Mar', 'varname_H06_I_Apr',
        'varname_H06_I_May', 'varname_H06_I_Jun', 'varname_H06_I_Jul', 'varname_H06_I_Aug',
        'varname_H06_I_Sep', 'varname_H06_I_Oct', 'varname_H06_I_Nov', 'varname_H06_I_Dec',
        'varname_H06_D_Jan', 'varname_H06_D_Feb', 'varname_H06_D_Mar', 'varname_H06_D_Apr',
        'varname_H06_D_May', 'varname_H06_D_Jun', 'varname_H06_D_Jul', 'varname_H06_D_Aug',
        'varname_H06_D_Sep', 'varname_H06_D_Oct', 'varname_H06_D_Nov', 'varname_H06_D_Dec',
        'varname_H06_purpose', 'varname_H06_I_mem_F', 'varname_H06_D_mem_F', 'varname_H06_I_mem_L',
        'varname_H06_D_mem_L', 'varname_hruContribIx', 'varname_hruContribId', 'varname_segId',
        'varname_segIndex', 'varname_downSegId', 'varname_downSegIndex', 'varname_upSegIds',
        'varname_upSegIndices', 'varname_rchOrder', 'varname_lakeId', 'varname_lakeIndex',
        'varname_isLakeInlet', 'varname_islake', 'varname_lakeModelType', 'varname_LakeTargVol',
        'varname_userTake', 'varname_goodBasin', 'varname_pfafCode', 'basRunoff', 'instRunoff',
        'dlayRunoff', 'sumUpstreamRunoff', 'KWTroutedRunoff', 'IRFroutedRunoff', 'KWroutedRunoff',
        'DWroutedRunoff', 'MCroutedRunoff', 'IRFvolume', 'KWTvolume', 'KWvolume', 'MCvolume',
        'DWvolume', 'KWfloodVolume', 'KWheight', 'MCfloodVolume', 'MCheight', 'DWfloodVolume',
        'DWheight', 'localSolute', 'soluteFlux', 'soluteMass', 'KWTinflow', 'IRFinflow',
        'KWinflow', 'MCinflow', 'DWinflow'
    }

   def readToml( self, infile, allowEmpty=False ):
       """
       Read and parse a mizuRoute TOML control file
       """
       logger.debug( "read in TOML file: "+infile )
       if ( not os.path.exists(infile) ):
          expect( False, "Input file to read does NOT exist: "+infile )

       with open( infile, "r" ) as ctlfile:
          content = ctlfile.read()

       parsed_toml = tomllib.loads(content)
       for key, val in parsed_toml.items():
          if key not in self.VALID_CONTROL_KEYS:
              expect( False, f"Unexpected variable in TOML control file: {key}" )
          val_str = str(val) if not isinstance(val, bool) else ('true' if val else 'false')
          self.set( key, val_str, allowNewName=True )

       if ( len(self.ctldict) == 0 and not allowEmpty ):
          expect( False, "No data was read from the file: "+infile )

       logger.debug( "File read" )
       self.fileRead = True

   def write_legacy( self, outfile ):
       """
       Write out a mizuRoute control file in legacy format
       """
       logger.debug( "Write out file: "+outfile )

       if ( os.path.exists(outfile) ):
          os.remove( outfile )
       ctlfile = open( outfile, "w" )
       vallen  = str(self.longestValue + 1)
       for name, value in self.ctldict.items():
          comment = self.comments.get(name, "")
          namelen = str(self.longestName - len(name) + 1)
          format = "<%s>%"+namelen+"s   %-"+vallen+"s    ! %s\n"
          ctlfile.write( format % (name, " ", value, comment) )

       ctlfile.close()

   def write( self, outfile ):
       """
       Write out a mizuRoute control file in TOML format
       """
       logger.debug( "Write out file: "+outfile )

       if ( os.path.exists(outfile) ):
          os.remove( outfile )
       ctlfile = open( outfile, "w" )
       for name, value in self.ctldict.items():
          val_str = str(value)
          if val_str.isdigit() or (val_str.startswith("-") and val_str[1:].isdigit()):
              formatted_val = val_str
          elif val_str.replace('.','',1).isdigit() or (val_str.startswith("-") and val_str[1:].replace('.','',1).isdigit()):
              formatted_val = val_str
          elif val_str.lower() in ['t', 'f', 'true', 'false', '.true.', '.false.']:
              formatted_val = 'true' if val_str.lower() in ['t', 'true', '.true.'] else 'false'
          else:
              formatted_val = f'"{val_str}"'

          comment = self.comments.get(name, "")
          comment_str = f" # {comment}" if comment else ""
          ctlfile.write( f"{name} = {formatted_val}{comment_str}\n" )

       ctlfile.close()

   def get( self, name ):
       """
       Return an element from the control file
       """
       return self.ctldict.get(name, "UNSET")

   def set( self, name, value, allowNewName=False, comment="" ):
       """
       Set an element in the control file
       """
       if ( len(name)  > self.longestName  ): self.longestName  = len(name)
       if ( len(str(value)) > self.longestValue ): self.longestValue = len(str(value))

       if ( not self._is_valid_name( name ) ):
          if ( allowNewName ):
             self.ctldict[name] = str(value)
             if comment:
                 self.comments[name] = comment
          else:
             expect( False, "set method is operating on a name that doesn't exist:"+name )
       else:
          self.ctldict[name] = str(value)
          if comment:
              self.comments[name] = comment

   def get_elmList( self ):
       """
       Get a copy of the list of elements in the file
       """
       if ( not self.is_read() ):
             expect( False, "mizuRoute control file was NOT read in yet, need to do that before returning list of elements" )

       return list(self.ctldict.keys())

   def _is_valid_name( self, name ):
       """
       Check if the name is valid
       """
       if ( self.is_read() ):
          return name in self.ctldict
       else:
          return False

   def is_read( self ):
       """
       Check if file has been read
       """
       return( self.fileRead )

#
# Unit testing for above classes
#
import unittest


class test_mizuRoute_control(unittest.TestCase):

   def setUp( self ):
       self.ctl = mizuRoute_control()

   def test_is_read( self ):
       self.assertFalse( self.ctl.is_read() )
       self.ctl.read( "SAMPLE_toml" )
       self.assertTrue( self.ctl.is_read() )

   def test_get_list_of_elments( self ):
       self.ctl.read( "SAMPLE_toml" )
       elist = self.ctl.get_elmList( )
       expected_subset = ['ancil_dir', 'input_dir', 'output_dir', 'sim_start', 'sim_end', 'fname_ntopOld',
                   'dname_sseg', 'dname_nhru',
                   'fname_ntopNew', 'seg_outlet', 'fname_qsim', 'vname_qsim',
                   'vname_time', 'vname_hruid', 'dname_xlon',
                   'dname_ylat', 'dname_time', 'dname_hruid', 'units_qsim', 'dt_qsim',
                   'is_remap', 'fname_remap', 'vname_hruid_in_remap',
                   'vname_weight', 'vname_qhruid', 'vname_num_qhru', 'dname_hru_remap',
                   'dname_data_remap', 'vname_i_index', 'vname_j_index',
                   'route_opt', 'is_flux_wm','fname_state_in',
                   'hydGeometryOption', 'topoNetworkOption',
                   'computeReachList', 'param_nml', 'varname_area', 'varname_length',
                   'varname_slope', 'varname_HRUid', 'varname_hruSegId',
                   'varname_segId', 'varname_downSegId']
       for expected_item in expected_subset:
           self.assertTrue(expected_item in elist, f"{expected_item} not in parsed list")

   def test_all_sample_fields_present( self ):
      ctl = mizuRoute_control.from_toml( "SAMPLE_toml" )
      self.assertTrue( ctl.is_read() )
      expected_set_keys = {
          'DWinflow', 'DWroutedRunoff', 'DWvolume', 'IRFinflow', 'IRFroutedRunoff', 'IRFvolume',
          'KWTinflow', 'KWTroutedRunoff', 'KWTvolume', 'KWinflow', 'KWroutedRunoff', 'KWvolume',
          'MCinflow', 'MCroutedRunoff', 'MCvolume', 'ancil_dir', 'basRunoff', 'case_name',
          'computeReachList', 'debug', 'desireId', 'dlayRunoff', 'dname_data_remap', 'dname_hru_remap',
          'dname_hruid', 'dname_nhru', 'dname_segid_wm', 'dname_sseg', 'dname_time', 'dname_time_wm',
          'dname_xlon', 'dname_ylat', 'doesBasinRoute', 'dt_qsim', 'dt_ro', 'fname_ntopNew',
          'fname_ntopOld', 'fname_qsim', 'fname_remap', 'fname_state_in', 'fname_wm',
          'histTimeStamp_offset', 'hw_drain_point', 'hydGeometryOption', 'input_dir', 'instRunoff',
          'is_flux_wm', 'is_lake_sim', 'is_remap', 'is_vol_wm', 'is_vol_wm_jumpstart',
          'newFileFrequency', 'ntopAugmentMode', 'outputFrequency', 'output_dir', 'param_nml',
          'restart_date', 'restart_write', 'route_opt', 'seg_outlet', 'sim_end', 'sim_start',
          'sumUpstreamRunoff', 'topoNetworkOption', 'units_qsim', 'varname_HRUid', 'varname_LakeTargVol',
          'varname_area', 'varname_downSegId', 'varname_hruSegId', 'varname_islake', 'varname_lakeModelType',
          'varname_length', 'varname_segId', 'varname_slope', 'vname_evapo', 'vname_flux_wm',
          'vname_hruid', 'vname_hruid_in_remap', 'vname_i_index', 'vname_j_index', 'vname_num_qhru',
          'vname_precip', 'vname_qhruid', 'vname_qsim', 'vname_segid_wm', 'vname_time',
          'vname_time_wm', 'vname_vol_wm', 'vname_weight'
      }
      parsed_keys = set(ctl.get_elmList())
      self.assertEqual( parsed_keys, expected_set_keys, f"Parsed keys do not match expected set. Diff: {parsed_keys ^ expected_set_keys}" )

   def test_unknown_toml_key_fails( self ):
      temp_file = "temp_unknown.toml"
      with open(temp_file, "w") as f:
          f.write('bogus_unknown_key = "invalid_value"\n')
      try:
          self.assertRaises( SystemExit, mizuRoute_control.from_toml, temp_file )
      finally:
          if os.path.exists(temp_file):
              os.remove(temp_file)

   def test_allow_empty( self ):
      self.ctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
      self.assertTrue( self.ctl.is_read() )

   def test_is_read_coupled( self ):
       self.assertFalse( self.ctl.is_read() )
       self.ctl.read( "SAMPLE-coupled_toml" )
       self.assertTrue( self.ctl.is_read() )

   def test_get_not_read( self ):
       value = self.ctl.get( "thing" )
       self.assertEqual( value, "UNSET" )

   def test_non_existant_file( self ):
       self.assertRaises( SystemExit, self.ctl.read, "file_does_NOT_EXIST.zztop" )

   def test_bad_file( self ):
       self.assertRaises( SystemExit, self.ctl.read, "README.md" )

   def test_get_after_set( self ):
       name = "thingwithlongname"
       value = "valuereturned"
       self.ctl.read( "SAMPLE_toml" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name )
       self.assertEqual( getvalue, value )

   def test_get_bad_name_after_set( self ):
       name = "thingwithlongname"
       name2 = name + "even_longer"
       value = "valuereturned"
       self.ctl.read( "SAMPLE_toml" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name2 )
       self.assertEqual( getvalue, "UNSET" )

   def test_set_doesnot_allow_newname( self ):
       name = "thingwithlongnamethatsnotonthefile"
       value = "valuetoset"
       self.ctl.read( "SAMPLE_toml" )
       self.assertRaises( SystemExit, self.ctl.set, name, value )

   def test_empty_file( self ):
       self.assertRaises( SystemExit, self.ctl.read, "../../cime_config/user_nl_mizuRoute" )

   def test_read_in_two_control_files( self ):
       self.ctl.read( "SAMPLE_toml" )
       newctl = mizuRoute_control()
       newctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
       self.assertEqual( [], newctl.get_elmList() )

   def test_write( self ):
       infile = "SAMPLE_toml"
       self.ctl.read( infile )
       outfile = "mizuRoute_in"
       self.ctl.write( outfile )
       self.assertTrue( os.path.exists(outfile) )
       os.remove( outfile )

   def test_read_legacy_control( self ):
       legacy_file = "temp_legacy_control"
       with open(legacy_file, "w") as f:
           f.write("<route_opt>        5   ! Legacy comment\n")
           f.write("<doesAccumRunoff>  1   ! Legacy comment\n")

       legacy_ctl = mizuRoute_control.from_control( legacy_file )
       self.assertEqual( legacy_ctl.get("route_opt"), "5" )
       self.assertEqual( legacy_ctl.get("doesAccumRunoff"), "1" )
       os.remove( legacy_file )

   def test_legacy_reader_fails_on_toml_data( self ):
       bad_control_file = "temp_toml_syntax_control"
       with open(bad_control_file, "w") as f:
           f.write("route_opt = 5\n")

       try:
           self.assertRaises( SystemExit, mizuRoute_control.from_control, bad_control_file )
       finally:
           if os.path.exists(bad_control_file):
               os.remove( bad_control_file )

   def test_toml_reader_fails_on_legacy_data( self ):
       bad_toml_file = "temp_legacy_syntax_toml"
       with open(bad_toml_file, "w") as f:
           f.write("<route_opt> 5 ! comment\n")

       try:
           self.assertRaises( (SystemExit, Exception), mizuRoute_control.from_toml, bad_toml_file )
       finally:
           if os.path.exists(bad_toml_file):
               os.remove( bad_toml_file )

   def test_factory_constructors( self ):
       toml_ctl = mizuRoute_control.from_toml( "SAMPLE_toml" )
       self.assertTrue( toml_ctl.is_read() )
       self.assertNotEqual( toml_ctl.get("route_opt"), "UNSET" )

       legacy_file = "temp_legacy_factory_control"
       with open(legacy_file, "w") as f:
           f.write("<route_opt>        3   ! Legacy comment\n")
       try:
           legacy_ctl = mizuRoute_control.from_control( legacy_file )
           self.assertEqual( legacy_ctl.get("route_opt"), "3" )
       finally:
           if os.path.exists(legacy_file):
               os.remove( legacy_file )

def main():
    import argparse
    parser = argparse.ArgumentParser(description="mizuRoute Control / TOML format conversion utility.")
    parser.add_argument("--from_toml", type=str, help="Input TOML control file path")
    parser.add_argument("--from_control", type=str, help="Input legacy control file path")
    parser.add_argument("--to_toml", type=str, help="Output TOML control file path")
    parser.add_argument("--to_control", type=str, help="Output legacy control file path")

    if len(sys.argv) == 1 or (len(sys.argv) > 1 and sys.argv[1].startswith("-v")):
        unittest.main()
        return

    args, unknown = parser.parse_known_args()

    if not args.from_toml and not args.from_control:
        unittest.main()
        return

    expect(not (args.from_toml and args.from_control), "Specify either --from_toml or --from_control, not both.")
    expect(not (args.to_toml and args.to_control), "Specify either --to_toml or --to_control, not both.")
    expect(args.to_toml or args.to_control, "Specify output target using either --to_toml=<fname> or --to_control=<fname>.")

    if args.from_toml:
        infile = args.from_toml
        ctl = mizuRoute_control.from_toml(infile)
    else:
        infile = args.from_control
        ctl = mizuRoute_control.from_control(infile)

    if args.to_toml:
        outfile = args.to_toml
        ctl.write(outfile)
    else:
        outfile = args.to_control
        ctl.write_legacy(outfile)

    print(f"Successfully converted {infile} -> {outfile}")

if __name__ == '__main__':
      main()
