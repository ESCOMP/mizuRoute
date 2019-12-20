"""
Python module to read in a sample mizuRoute control file into a data type
as an array of keys as well as a dictionary of the values. The Dictionary
can then be modified and output as a new file.

Erik Kluzek
"""

import sys

sys.path.append( "../../cime/scripts/lib" );

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type

import six

logger = logging.getLogger(__name__)

class mizuRoute_control(object):

   """ Object to hold a dictionary of settings for mizuRoute control """

   # Class Data:
   fileread = False

   def read( self, infile ):
       """
       Read and parse a mizuRoute control file
       """
       self.fileread = True


   def write( self, outfile ):
       """
       Write out a mizuRoute control file
       """

   def get( self, name ):
       """
       Return an element from the control file
       """

   def set( self, name ):
       """
       Set an element in the control file
       """


   def __is_valid_name( self, name ):
       """
       Check if the name is valid
       """
       if ( self.is_read() ):
          return( True )
       else:
          return( False )

   def is_read( self ):
       """
       Check if file has been read
       """
       return( self.fileread )

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

if __name__ == '__main__':
     unittest.main()
