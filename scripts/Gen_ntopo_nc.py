#!/usr/bin/env python

"""
 Generate Network Topology netCDF from hru shapefile and segment shapefile
 Required input data for this script.
 1. hru shapefile (polygon)
 2. stream segment shapefile (line)
 3. shapefile meta data 


Usage: 
   ./Gen_ntopo_nc.py -[h] -c<in_hru_shp> -s<in_seg_shp> -a<att_meta> -o<outnc>
    -c  input hru catchment shapefle 
    -s  input segment shapefile 
    -a  metadata for attributes 
    -o  output netcdf
    -h  help

"""

# Import modules
import sys 
import time
import os
import getopt
import gdal,ogr,osr 
import numpy as np
import math 
import netCDF4 as nc4

#-----------------------------------------------
#  HARD CODED VARIABLES
#-----------------------------------------------
GISDRV  = ogr.GetDriverByName('ESRI Shapefile')
SEGIDNM = 'seg_id'
HRUIDNM = 'hruid'
#-----------------------------------------------
# Classes
#-----------------------------------------------
class NTOPO: 
    """ topology object - combination of segment and hru """
    def __init__(self, inseg, inhru, inatm, ncout):
      #Input
      self.inseg = inseg  # stream segment shapefile name
      self.inhru = inhru  # hru shapefile name
      self.inatm = inatm  # shapefile meta data name
      self.ncout = ncout  # netcdf output name

      # Initialize list and dictionary
      self.segID = []
      self.hruID = []
      self.meta   = AutoVivification() 
      self.attVal = AutoVivification() 

    def process(self):
      """Process everything """
      self.getAttMeta()  
      self.getSegID()  
      self.getHruID()  
      self.getAtt()
      self.compAtt()
      self.write_nc()

    def getAttMeta(self):
      """ read metadata for hru and segment shp attributes """
      with open(self.inatm) as f:
        cnt=0
        while True :
          cols=f.readline().split(',')
          if len(cols)<2:
            break
          self.meta[cols[2]]['shptyp']    = cols[0].strip()
          self.meta[cols[2]]['how2get']   = cols[1].strip()
          self.meta[cols[2]]['dtyp']      = cols[3].strip()
          self.meta[cols[2]]['long_name'] = cols[4].strip()
          self.meta[cols[2]]['units']     = cols[5].strip()

    def getSegID(self):
      """ read segment ID from segment shapefile"""
      self.segID=getShpID(self.inseg, SEGIDNM) 
      
    def getHruID(self):
      """ read hru ID from hru shapefile"""
      self.hruID=getShpID(self.inhru, HRUIDNM) 

    def getAtt(self):
      """ read segment or hru attributes into attVal nested dictionary """
      # nested dictionary -  attVal[attribute name] {id;value}
      for attn, val in self.meta.iteritems():
        if val['shptyp'] == 'seg' and val['how2get'] == 'read':
          tmp = getShpField(self.inseg,SEGIDNM,attn)  # tmp  {shpid: value} for current attribute 
          for atkey, atval in tmp.iteritems():
            self.attVal[attn][atkey] = atval 
        if val['shptyp'] == 'hru' and val['how2get'] == 'read':
          tmp = getShpField(self.inhru,HRUIDNM,attn) 
          for atkey, atval in tmp.iteritems():
            self.attVal[attn][atkey] = atval 

    def compAtt(self):
      """ comp segment attributes AttVal dictionary """
      # attribute name is now hard coded so make sure these are the same as attribute meta file
      for attn, val in self.meta.iteritems():
        if attn == 'Drop':
          for sids in self.segID:
            self.attVal[attn][sids] = self.attVal['TopElev'][sids]-self.attVal['BotElev'][sids]
        if attn == 'Slope':
          for sids in self.segID:
            self.attVal[attn][sids] = tmpval = (self.attVal['TopElev'][sids]-self.attVal['BotElev'][sids])/self.attVal['Length'][sids]

    def write_nc(self):
      """ Creates an netcdf file containing topology info """ 
      # Create netcdf file for this simulation
      rootgrp = nc4.Dataset(self.ncout, 'w', format='NETCDF4')
      # Dimension size 
      hruSize=getShpSize(self.inhru)
      segSize=getShpSize(self.inseg)
      # Create dimensions
      seg_id = rootgrp.createDimension(SEGIDNM, segSize)           # Length = number of seg_id values
      hru_id = rootgrp.createDimension(HRUIDNM, hruSize)           # Length = number of hru_id values
      # Create fixed-length variables
      seg_ids = rootgrp.createVariable(SEGIDNM,'i4',(SEGIDNM))     # Coordinate Variable (32-bit signed integer)
      hru_ids = rootgrp.createVariable(HRUIDNM,'i4',(HRUIDNM))     # Coordinate Variable (32-bit signed integer)
      # Fill in coordinate variables 
      seg_ids[:]        = self.segID
      hru_ids[:]        = self.hruID
      # go through segment
      for attn, val in self.meta.iteritems(): 
        if val['shptyp'] == 'seg':
          ncvar           = rootgrp.createVariable(attn, self.meta[attn]['dtyp'], (SEGIDNM))
          #Make sure order of attribute value is in the order of id
          sval=[]
          for sid in self.segID:
            sval.append(self.attVal[attn][sid])
          ncvar[:]        = sval 
        elif val['shptyp'] == 'hru':
          ncvar           = rootgrp.createVariable(attn, self.meta[attn]['dtyp'], (HRUIDNM))
          hval=[]
          for hid in self.hruID:
            hval.append(self.attVal[attn][hid])
          ncvar[:]        = hval 
        #attributes for variable
        ncvar.long_name = self.meta[attn]['long_name']  # long name
        ncvar.units     = self.meta[attn]['units']      # units
      # Fill in global attributes
      rootgrp.history = 'Created on %s' %time.ctime()
      rootgrp.history = 'Created with %s' %__file__
      # Close file
      rootgrp.close()

class AutoVivification(dict):
  """Implementation of perl's autovivification feature to initialize structure."""
  def __getitem__(self, item):
    try:
      return dict.__getitem__(self, item)
    except KeyError:
      value = self[item] = type(self)()
    return value

#-----------------------------------------------
# Modules
#-----------------------------------------------
def getShpSize(inShp):
    """ Get number of features in shapefile """
    shp = GISDRV.Open(inShp,0)
    layer = shp.GetLayer()
    
    nFtres = layer.GetFeatureCount()

    return nFtres 

def getShpID(inShp, idnm):
    """ Get shapefile ID in list """
    # Open the shapefile file with GDAL, with read only access
    shp = GISDRV.Open(inShp,0)
    layer = shp.GetLayer()
    
    nFtres = layer.GetFeatureCount()

    ID=[]
    for i in range(nFtres):
      ftre = layer.GetFeature(i)
      ID.append(ftre.GetField(idnm))

    return ID 

def getShpField(inShp, shp_id, fldname):
    """ Get shapefile field value with dictionary data format hru_field {shp_id; fieldvalue} """
    # Open the shapefile file with GDAL, with read only access
    shp = GISDRV.Open(inShp,0)
    layer = shp.GetLayer()
    
    nFtres = layer.GetFeatureCount()

    fldval={}
    for i in range(nFtres):
      ftre = layer.GetFeature(i)
      sid = ftre.GetField(shp_id)
      fldval[sid] = ftre.GetField(fldname)

    return fldval

##############################################################
use = '''
Usage: %s -[h] -c<in_hru_shp> -s<in_seg_shp> -a<att_meta> -o<outnc>

    -c  input hru catchment shapefle 
    -s  input segment shapefile 
    -a  metadata for attributes 
    -o  output netcdf
    -h  help

'''
if __name__ == '__main__':

  def usage():
    sys.stderr.write(use % sys.argv[0])
    sys.exit(1)

  try:
    (opts, args) = getopt.getopt(sys.argv[1:], 'hc:s:a:o:')
  except getopt.error:
    usage()

  inhru = None
  inseg = None
  inatm = None
  for (opt,val) in opts:
    if opt == '-h':
      usage()
    elif opt == '-c':
      inhru = val
    elif opt == '-s':
      inseg = val
    elif opt == '-a':
      inatm = val
    elif opt == '-o':
      outnc = val
    else:
      usage()

  if not inhru or not inseg or not inatm:
    usage()

  # Gather the size of each dataset
  a=NTOPO(inseg, inhru, inatm, outnc) 
  a.process()

  print "completed processing...." 
