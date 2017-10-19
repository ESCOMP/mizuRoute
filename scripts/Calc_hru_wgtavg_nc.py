#!/usr/bin/env python

#this is NOT meant to be anything more than a start for people.
#it is NOT a general use utility.
#no one involved in creating or distributing this script has any
#responsibility for its function or output.

import sys
import os
import time 
import getopt
import numpy as np
import math 
import netCDF4 as nc4

############################################ 
#         hardcoded variables              #
############################################
# Name of netCDF variable for polygon ID in input data
IDNM='hru_id'               
# Name of netCDF variable for weight in weight data
WGTNM='weight'
OVRPLYNM='overlapPolyId'
LATNM='latitude'
LONNM='longitude'
OVRNM='overlaps'
############################################ 
#              Class                       #
############################################
class wgtnc:
    """object of basin netCDF including hru and areal weight/lat/lon of others polygons  """
    def __init__(self,ncName):
        """Initialization """
        self.ncName=ncName

    def getWgtHru(self,hru):
        """For given hru id, get weight of the intersected polygons and associated ID and lat/lon"""
        wgtAll        = getNetCDFData(self.ncName, WGTNM) 
        overlapsIdAll = getNetCDFData(self.ncName, OVRPLYNM) 
        latAll        = getNetCDFData(self.ncName, LATNM)
        lonAll        = getNetCDFData(self.ncName, LONNM)
        overlapsAll   = getNetCDFData(self.ncName, OVRNM)

        self.hruList = self.getHruID()            # get hru id list
        idx=self.hruList.index(hru)               # get indix in array corresponding hru
        self.wgt        = list(wgtAll[idx])       # Get overlapping poly's wgt list for hru 
        self.overlapsId = list(overlapsIdAll[idx])# Get overlapping poly's wgt list for hru 
        self.lat        = list(latAll[idx])       # Get overlapping poly's lat list for hru
        self.lon        = list(lonAll[idx])       # Get overlapping poly's lon list for hru
        self.overlaps   = overlapsAll[idx]        # Get number of overlapping polys
        
        return (self.wgt, self.overlapsId, self.lat, self.lon, self.overlaps)

    def getHruIdName(self):
        """ get Name of hru ID 
            it must be the 1st dimension """
        f = nc4.Dataset(self.ncName,'r')
        self.dim = f.dimensions
        self.dimName = self.dim.keys()
        self.hruIdName = self.dimName[0]

        return self.hruIdName
        
    def getHruID(self):
        """ get hru ID list of basin"""
        self.HruIdName=self.getHruIdName()
        self.hruid = list(getNetCDFData(self.ncName, self.HruIdName))
        return self.hruid

############################################
#            Modules                       #
############################################
def getNetCDFData(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
    return data

def getNetCDFAtt(fn, varname,attName):
    """Read attribute of <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    var = f.variables[varname]
    attData = getattr(var,attName)
    return attData

def writeNetCDFData(fn, var, varname):
  """ Write <var>[time,hru] array in netCDF4 file,<fn> and variabl of <varname> """
  ncfile = nc4.Dataset(fn,'w',format='NETCDF4')

  # Data should be 2D [time x hru]
  dim1size=var.shape[0]
  dim2size=var.shape[1]
 
  dim_1 = ncfile.createDimension('dim1',dim1size )  # hru axis
  dim_2 = ncfile.createDimension('dim2',dim2size )  # hru axis
#  dim_1 = ncfile.createDimension('time',None )  # time axis - record dimension

  # Define a 2D variable to hold the var
  val = ncfile.createVariable(varname,'f4',('dim1','dim2'))
  
  # Write grid 
  val[:,:]=var

  # Write basic global attribute
  ncfile.history = 'Created ' + time.ctime(time.time())
  ncfile.source = os.path.dirname(os.path.abspath(__file__))+__file__[1:]
  
  ncfile.close()

def compAvgVal(nc_wgt,nc_in,varname):
  """Compute areal weighted avg value of <varname> in <nc_in> for each hru based on hru's weight in <nc_wgt>""" 
  wgt = wgtnc(nc_wgt)             #instantaneous of wgtnc object
  hruIdName = wgt.getHruIdName()  # Get precise name of variable for hru id 
  hruIDs    = wgt.getHruID()      #get hruID list 

  dataVal  = getNetCDFData(nc_in,varname) # Get data value 
  FillVal  = getNetCDFAtt(nc_in,varname,'_FillValue') # Get data value 
  IdVal    = getNetCDFData(nc_in,IDNM) 
  dim1size = dataVal.shape[0]

  #Initialize wgtVal[ntime,nhru] 
  wgtVal = np.zeros((dim1size,len(hruIDs)))
  #Loop through each hru polygon 
  for i in range(len(hruIDs)):
    print "Computing weighted values over hru%d" %hruIDs[i]
    # Get list of wgt, lat, and lon for corresponding hru
    (wgtval, overlapsId, lat, lon, overlaps)=wgt.getWgtHru(hruIDs[i])
    wgtArray = np.asarray(wgtval)
    
    #Count missing cells (No values in original input polygon) and valid cells
    numvoid = 0
    # Go through wgt list and replace value with zero for following cases
    # where overlapping polygon has missing value - case1
    # where overlapping polygon is outside nc_wgt domain - case2
    for j in range(len(wgtArray)): # Go through each overlapping polygon

      if overlapsId[j] >= 0:              # if there is at least one overlapping polygon
        # find index of grid cell that match up with hru id of overlapsId[j] 
        ij=np.where(IdVal==overlapsId[j])
        # if nc_in netCDF does not cover nc_wgt domain - case2  
        if len(ij[0]) == 0:
          numvoid = numvoid+1
          wgtArray[j]=0
        else:
          a = dataVal[10,ij[0][0],ij[1][0]].tolist()
          #if value of overlapping polygon is missing data -case1
          if a == None: 
            numvoid = numvoid+1
            wgtArray[j] = 0 
    
    Val2 = np.zeros((dim1size))
    if np.nansum(wgtArray) > 0.0:
      # Adjust weight value if valid weight value (> 0) exist in list
      newWgtArray = [x/np.nansum(wgtArray) for x in wgtArray]
      # Go through each overlapping polygon
      c=0
      for j in range(len(newWgtArray)): 
        if overlapsId[j] >= 0: # if there is at least one intersecting polygon

          ij=np.where(IdVal==overlapsId[j])  
          yj = ij[0]
          xi = ij[1]
          Val1=newWgtArray[j]*dataVal[:,yj,xi]
          #just check If selected hru id from runoff nc matches w/ wieght nc, if not,something wrong exit 
          if (IdVal[yj,xi] != overlapsId[j]):
             print "hru_id %d not exist, Sayonara" %IdVal[yj,xi]
             sys.exit()
          # Concatenate time series array, Val1 to oriVal containin time serise from all the intersected polygon 
          if Val1.size: #if Val1 is not empty i.e., exclude empty value in weight list
            c += 1
            if c == 1:
              oriVal = Val1
            else:
              oriVal=np.concatenate((oriVal, Val1), axis=1)
      # This is weighted average value at current hru
      Val2 = oriVal.sum(axis=1)
      del oriVal

    #if there is no valid weight values (i.e. all zero)- nothing overlapping for current hru, so missing value
    else:
      Val2 = np.ones((dim1size))*-9999.0
    # store weighted average time series, Val2 in final array wgtVal
    wgtVal[:,i]= Val2
    
    del Val2

  return wgtVal

############################################
#                Main                      #
############################################
use = '''
Usage: %s -[h] <weight_netCDF> <input_netCDF> <variable_name_in_input_netCDF> <output_netCDF>
        -h  help
'''
if __name__ == '__main__':

    def usage():
        sys.stderr.write(use % sys.argv[0])
        sys.exit(1)
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], 'h')
    except getopt.error:
        usage()

    verbose = False
    grid_info = False
    proj_info = True
    for (opt,val) in opts:
        if opt == '-h':
            usage()
        elif opt == '-v':
            verbose = True
        else:
            raise OptionError, opt
            usage()

    if len(args) == 4:
      # Read three argument
      nc_wgt  = args[0]
      nc_in   = args[1]  
      varname = args[2]  
      nc_out  = args[3]

      wgtVal=compAvgVal(nc_wgt,nc_in,varname)
      writeNetCDFData(nc_out, wgtVal, varname)     
    else:
      usage()

