#!/usr/bin/env python

#this is NOT meant to be anything more than a start for people.
#it is NOT a general use utility.
#no one involved in creating or distributing this script has any
#responsibility for its function or output.

import sys
import os
import multiprocessing as mp
from functools import partial
import getopt
import xarray as xr
import numpy as np
import netCDF4 as nc4

############################################ 
#         hardcoded variables              #
############################################
# number of processors (= number of chunks of hru to be processed)
NPROC=4
# Name of netCDF variable for polygon ID in input data
HDIM='hru_id'
TDIM='Time'
# Name of netCDF variable for weight in weight data
WGTNM='weight'
OVRPLYNM='overlapPolyId'
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
        overlapsAll   = getNetCDFData(self.ncName, OVRNM)

        self.hruList = self.getHruID()            # get hru id list
        idx=self.hruList.index(hru)               # get indix in array corresponding hru
        self.wgt        = np.asarray(list(wgtAll[idx]))       # Get overlapping poly's wgt list for hru
        self.overlapsId = np.asarray(list(overlapsIdAll[idx])) # Get overlapping poly's wgt list for hru
        self.overlaps   = overlapsAll[idx]        # Get number of overlapping polys
        return (self.wgt, self.overlapsId, self.overlaps)

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

def compAvgVal(nc_wgt, nc_in, varname, chunk):
  """Compute areal weighted avg value of <varname> in <nc_in> for each hru based on hru's weight in <nc_wgt>""" 
  wgt       = wgtnc(nc_wgt)       #instantaneous of wgtnc object
  hruIDs    = wgt.getHruID()      #get hruID list 
  hrudata   = split(hruIDs,NPROC)
  hrulist   = hrudata[chunk]

  dataVal  = getNetCDFData(nc_in,varname) # Get data value 
  FillVal  = getNetCDFAtt(nc_in,varname,'_FillValue') # Get data value 
  ghruid   = getNetCDFData(nc_in,HDIM)
  dim1size = dataVal.shape[0]
  #Initialize wgtVal[ntime,nhru] 
  wgtVal = np.ones((dim1size,len(hrulist)))*-9999.0
  #Loop through each hru polygon 
  for i in range(len(hrulist)):
    print hrulist[i]
    # Get list of wgt for corresponding hru
    (wgtArray, overlapsId, overlaps)=wgt.getWgtHru(hrulist[i])
    #Count missing cells (No values in original input polygon) and valid cells
    numvoid = 0
    # Go through wgt list and replace value with zero for following cases
    # where overlapping polygon has missing value - case1
    # where overlapping polygon is outside nc_wgt domain - case2
    a = np.zeros( (dim1size,len(wgtArray)) )
    for j in range(overlaps):                     # Go through each overlapping polygon
      # find index of grid cell that match up with hru id of overlapsId[j]
      ij=np.column_stack(np.where(ghruid==overlapsId[j]))
      # if nc_in netCDF does not cover nc_wgt domain - case2
      if len(ij[0]) == 0:
        numvoid+=1
        wgtArray[j]=0.0
      else:
        a[:,j] = dataVal[:,ij[0,0],ij[0,1]]
        #if value of overlapping polygon is missing data -case1
        if any(a[:,j]==FillVal):
          numvoid+=1
          wgtArray[j] = 0.0
    if np.nansum(wgtArray) > 0.0:
      # Adjust weight value if valid weight value (> 0) exist in list
      wgtArray = wgtArray/np.nansum(wgtArray)
      wgtArray = np.where(np.isnan(wgtArray),0,wgtArray)
      wgtVal[:,i]=np.dot(a,wgtArray.T).T

  return wgtVal

def split(hrulist,n):
  subsize=len(hrulist)/n
  idx1=0
  hrudata=[]
  for i in range(n-1):
    idx2=idx1+subsize
    hrudata.append(hrulist[idx1:idx2])
    idx1=idx2
  hrudata.append(hrulist[idx1:])
  return hrudata

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

      func=partial(compAvgVal, nc_wgt, nc_in, varname)

      pool = mp.Pool(processes=NPROC,maxtasksperchild=1)
      results=pool.map(func,range(NPROC))
      result=np.column_stack(results)

      with xr.open_dataset(nc_in) as ds:
        timedata  = ds.coords[TDIM]     # Get time coordinate variable
      with xr.open_dataset(nc_wgt) as ds:
        hrudata = ds[HDIM].values

      encoding={varname: {'dtype': 'float32', '_FillValue': -9999.0}}
      foo = xr.DataArray(result, coords=[('time',timedata), ('hru',hrudata)], name=varname)
      foo.encoding = encoding[varname]
      foo.to_netcdf('test.nc')
    else:
      usage()

