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
   fileRead = False
   dict = {}
   keyList = []

   def read( self, infile ):
       """
       Read and parse a mizuRoute control file
       """
       self.fileRead = True


   def write( self, outfile ):
       """
       Write out a mizuRoute control file
       """

   def get( self, name ):
       """
       Return an element from the control file
       """
       if ( self.__is_valid_name( name ) ):
          return( self.dict[name] )
       else:
          return( "UNSET" )

   def set( self, name, value ):
       """
       Set an element in the control file
       """
       if ( not self.__is_valid_name( name ) ):
          self.keyList.append(name)

       self.dict[name] = value


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

   def test_get_not_read( self ):
       value = self.ctl.get( "thing" )
       self.assertEqual( value, "UNSET" )

   def test_get_after_set( self ):
       name = "thingwithlongname"
       value = "valuereturned"
       self.ctl.read( "SAMPLE.control" )
       self.ctl.set( name, value )
       getvalue = self.ctl.get( name )
       self.assertEqual( getvalue, value )

   def test_get_bad_name_after_set( self ):
       name = "thingwithlongname"
       name2 = name + "even_longer"
       value = "valuereturned"
       self.ctl.read( "SAMPLE.control" )
       self.ctl.set( name, value )
       getvalue = self.ctl.get( name2 )
       self.assertEqual( getvalue, "UNSET" )


if __name__ == '__main__':
     unittest.main()
