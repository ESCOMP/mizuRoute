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
   lineMatch = '^<(.+?)>\s+(\S+)\s+\!(.+)$' # Pattern to match for legacy lines
   longestName = 0                          # Longest name
   longestValue = 0                         # Longest value

   def __init__(self):
      self.ctldict = collections.OrderedDict()     # Ordered dictionary of control elements
      self.comments = {}                           # Comments associated with keys

   def read( self, infile, allowEmpty=False ):
       """
       Read and parse a mizuRoute control file
       """
       if ( infile.endswith(".toml") ):
           return self.readToml( infile, allowEmpty=allowEmpty )

       logger.debug( "read in file: "+infile )
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
       self.ctl.read( "SAMPLE.toml" )
       self.assertTrue( self.ctl.is_read() )

   def test_get_list_of_elments( self ):
       self.ctl.read( "SAMPLE.toml" )
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

   def test_allow_empty( self ):
       self.ctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
       self.assertTrue( self.ctl.is_read() )

   def test_is_read_coupled( self ):
       self.assertFalse( self.ctl.is_read() )
       self.ctl.read( "SAMPLE-coupled.toml" )
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
       self.ctl.read( "SAMPLE.toml" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name )
       self.assertEqual( getvalue, value )

   def test_get_bad_name_after_set( self ):
       name = "thingwithlongname"
       name2 = name + "even_longer"
       value = "valuereturned"
       self.ctl.read( "SAMPLE.toml" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name2 )
       self.assertEqual( getvalue, "UNSET" )

   def test_set_doesnot_allow_newname( self ):
       name = "thingwithlongnamethatsnotonthefile"
       value = "valuetoset"
       self.ctl.read( "SAMPLE.toml" )
       self.assertRaises( SystemExit, self.ctl.set, name, value )

   def test_empty_file( self ):
       self.assertRaises( SystemExit, self.ctl.read, "../../cime_config/user_nl_mizuRoute" )

   def test_read_in_two_control_files( self ):
       self.ctl.read( "SAMPLE.toml" )
       newctl = mizuRoute_control()
       newctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
       self.assertEqual( [], newctl.get_elmList() )

   def test_write( self ):
       infile = "SAMPLE.toml"
       self.ctl.read( infile )
       outfile = "mizuRoute_in"
       self.ctl.write( outfile )
       self.assertTrue( os.path.exists(outfile) )
       os.remove( outfile )

   def test_read_legacy_control( self ):
       legacy_file = "temp_legacy.control"
       with open(legacy_file, "w") as f:
           f.write("<route_opt>        5   ! Legacy comment\n")
           f.write("<doesAccumRunoff>  1   ! Legacy comment\n")

       legacy_ctl = mizuRoute_control()
       legacy_ctl.read( legacy_file )
       self.assertEqual( legacy_ctl.get("route_opt"), "5" )
       self.assertEqual( legacy_ctl.get("doesAccumRunoff"), "1" )
       os.remove( legacy_file )

if __name__ == '__main__':
     unittest.main()
