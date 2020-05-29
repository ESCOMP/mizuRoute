"""
Python module to read in a sample mizuRoute control file into a data type
as an array of keys as well as a dictionary of the values. The Dictionary
can then be modified and output as a new file.

Erik Kluzek
"""

import sys, re

sys.path.append( "../../cime/scripts/lib" );
sys.path.append( "../../../../cime/scripts/lib" );

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type, run_cmd_no_fail

import six

logger = logging.getLogger(__name__)

class mizuRoute_control(object):

   """ Object to hold a dictionary of settings for mizuRoute control """

   # Class Data:
   fileRead = False                         # If file has been read or not
   lineMatch = '^<(.+?)>\s+(\S+)\s+\!(.+)$' # Pattern to match for lines
   longestName = 0                          # Longest name
   longestValue = 0                         # Longest value

   def __init__(self):
      self.ctldict = {}                             # Dictionary of control elments
      self.keyList = []                             # List of keys for control elements
      self.lines = []                               # Lines of the entire file read in

   def read( self, infile, allowEmpty=False ):
       """
       Read and parse a mizuRoute control file
       """
       # Read the whole file and save each line as object data
       logger.debug( "read in file: "+infile )
       if ( not os.path.exists(infile) ):
          expect( False, "Input file to read does NOT exist: "+infile )

       ctlfile = open( infile, "r" )
       self.lines = ctlfile.readlines()
       ctlfile.close()

       # Loop through each line in the file
       for line in self.lines:
          # Ignore comment lines
          if ( not line.find( "!" ) == 0 ):
             match = re.search( self.lineMatch, line )
             if ( not match ):
                expect( False, "Error in reading in line:"+line )
             else:
                name = match.group(1)
                value = match.group(2)
                self.set( name, value, allowNewName=True )


       # If no data was read -- abort with an error
       if ( len(self.keyList) == 0 and not allowEmpty ):
          expect( False, "No data was read from the file: "+infile )

       # Mark the file as read
       logger.debug( "File read" )
       self.fileRead = True


   def write( self, outfile ):
       """
       Write out a mizuRoute control file
       """
       logger.debug( "Write out file: "+outfile )

       if ( os.path.exists(outfile) ):
          os.remove( outfile )
       ctlfile = open( outfile, "w" )
       vallen  = str(self.longestValue + 1)
       # Loop through each line in the file
       for line in self.lines:
          # Write comment lines as is
          if ( line.find( "!" ) == 0 ):
             ctlfile.write( line )
          else:
             match = re.search( self.lineMatch, line )
             if ( not match ):
                expect( False, "Error in for output line:"+line )
             name = match.group(1)
             value = self.get( name )
             comment = match.group(3)
             namelen = str(self.longestName - len(name) + 1)
             format = "<%s>%"+namelen+"s   %-"+vallen+"s    ! %s\n"
             ctlfile.write( format % (name, " ", value, comment) )

       ctlfile.close()

   def get( self, name ):
       """
       Return an element from the control file
       """
       if ( self.__is_valid_name( name ) ):
          return( self.ctldict[name] )
       else:
          return( "UNSET" )

   def set( self, name, value, allowNewName=False ):
       """
       Set an element in the control file
       """
       self.ctldict[name] = value
       # Check for longest value and name
       if ( len(name)  > self.longestName  ): self.longestName  = len(name)
       if ( len(value) > self.longestValue ): self.longestValue = len(value)

       if ( not self.__is_valid_name( name ) ):
          if ( allowNewName ):
             self.keyList.append(name)
          else:
             expect( False, "set method is operating on a name that doesn't exist:"+name )

   def get_elmList( self ):
       """
       Get a copy of the list of elements in the file
       """
       if ( not self.is_read() ):
             expect( False, "mizuRoute control file was NOT read in yet, need to do that before returning list of elements" )

       elmList = list(self.keyList)
       return( elmList )

   def __is_valid_name( self, name ):
       """
       Check if the name is valid
       """
       if ( self.is_read() ):
          try:
             idx =  self.keyList.index(name)
             return( True )
          except  ValueError:
             return( False )
       else:
          return( False )

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
       self.ctl.read( "SAMPLE.control" )
       self.assertTrue( self.ctl.is_read() )

   def test_get_list_of_elments( self ):
       self.ctl.read( "SAMPLE.control" )
       elist = self.ctl.get_elmList( )
       expected = ['ancil_dir', 'input_dir', 'output_dir', 'sim_start', 'sim_end', 'fname_ntopOld',
                   'dname_sseg', 'dname_nhru',
                   'fname_ntopNew', 'seg_outlet', 'fname_qsim', 'vname_qsim',
                   'vname_time', 'vname_hruid', 'dname_xlon',
                   'dname_ylat', 'dname_time', 'dname_hruid', 'units_qsim', 'dt_qsim',
                   'is_remap', 'fname_remap', 'vname_hruid_in_remap',
                   'vname_weight', 'vname_qhruid', 'vname_num_qhru', 'dname_hru_remap',
                   'dname_data_remap', 'vname_i_index', 'vname_j_index',
                   'route_opt', 'fname_state_in',
                   'hydGeometryOption', 'topoNetworkOption',
                   'computeReachList', 'param_nml', 'varname_area', 'varname_length',
                   'varname_slope', 'varname_HRUid', 'varname_hruSegId',
                   'varname_segId', 'varname_downSegId']
       self.assertEqual( expected, elist )

   def test_allow_empty( self ):
       self.ctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
       self.assertTrue( self.ctl.is_read() )

   def test_is_read_coupled( self ):
       self.assertFalse( self.ctl.is_read() )
       self.ctl.read( "SAMPLE-coupled.control" )
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
       self.ctl.read( "SAMPLE.control" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name )
       self.assertEqual( getvalue, value )

   def test_get_bad_name_after_set( self ):
       name = "thingwithlongname"
       name2 = name + "even_longer"
       value = "valuereturned"
       self.ctl.read( "SAMPLE.control" )
       self.ctl.set( name, value, allowNewName=True )
       getvalue = self.ctl.get( name2 )
       self.assertEqual( getvalue, "UNSET" )

   def test_set_doesnot_allow_newname( self ):
       name = "thingwithlongnamethatsnotonthefile"
       value = "valuetoset"
       self.ctl.read( "SAMPLE.control" )
       self.assertRaises( SystemExit, self.ctl.set, name, value )


   def test_empty_file( self ):
       self.assertRaises( SystemExit, self.ctl.read, "../../cime_config/user_nl_mizuRoute" )

   def test_read_in_two_control_files( self ):
       # Read in two control files make sure their list of elements is different
       self.ctl.read( "SAMPLE.control" )
       newctl = mizuRoute_control()
       newctl.read( "../../cime_config/user_nl_mizuRoute", allowEmpty=True )
       self.assertEqual( [], newctl.get_elmList() )

   def test_write( self ):
       infile = "SAMPLE.control"
       self.ctl.read( infile )
       outfile = "mizuRoute_in"
       self.ctl.write( outfile )
       if ( not run_cmd_no_fail( "diff -wb "+infile+" "+outfile ) == "" ):
          expect( False, "Write of input file results in something different" )
       os.remove( outfile )

if __name__ == '__main__':
     unittest.main()
