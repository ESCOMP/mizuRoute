#!/usr/bin/env python

#this is NOT meant to be anything more than a start for people.
#it is NOT a general use utility.
#no one involved in creating or distributing this script has any
#responsibility for its function or output.

# Created 
# 1/15/14 
# Updated 
# 2/25/14  

import sys
import os
import time 
import getopt
import numpy as np
import math 
import netCDF4 as nc4

############################################ 
#              Class                       #
############################################
class wgtnc:
    """object of basin netCDF including hru and areal weight/lat/lon of others polygons  """
    def __init__(self,ncName):
        """Initialization """
        self.ncName=ncName

    def getWgtHru(self,hru):
        """For given hru id, get weight of the intersected polygons and associated lat/lon"""
        wgtAll        = getNetCDFData(self.ncName, 'weight') 
        overlapsIdAll = getNetCDFData(self.ncName, 'overlapPolyId') 
        latAll        = getNetCDFData(self.ncName, 'latitude')
        lonAll        = getNetCDFData(self.ncName, 'longitude')
        overlapsAll   = getNetCDFData(self.ncName, 'overlaps')

        self.hruList = self.getHruID() # get hru id list
        Index=self.hruList.index(hru)  # get indix in array corresponding hru
        self.wgt        = list(wgtAll[Index]) #Get overlapping poly's wgt list for hru 
        self.overlapsId = list(overlapsIdAll[Index]) #Get overlapping poly's wgt list for hru 
        self.lat        = list(latAll[Index]) #Get overlapping poly's lat list for hru
        self.lon        = list(lonAll[Index]) #Get overlapping poly's lon list for hru
        self.overlaps   = overlapsAll[Index] #Get number of overlapping polys
        
        return (self.wgt, self.overlapsId, self.lat, self.lon, self.overlaps)

    def getHruIdName(self):
        """ get Name of hru ID """
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
#  dim_2 = ncfile.createDimension('hru',dim2size )  # hru axis

#  # Define time dimension variable to hold the var
#  t   = ncfile.createVariable('time','f4',('time'))

  # Define a 2D variable to hold the var
  val = ncfile.createVariable(varname,'f4',('dim1','dim2'))
#  val = ncfile.createVariable(varname,'f4',('time','hru'))
  
  # Write grid 
#  t[:]=times   # NM added
  val[:,:]=var

#  ### NM added
#  #Attribute    
#  # time axis
#  t.units = 'days since 1950-01-01 00:00:00'
#  t.calendar = 'gregorian'
#  t.long_name = 'Time axis'
#  t.time_origin = '1950-JAN-01 00:00:00'
#  # hru axis
#  val.units = 'mm/month'
#  val.associate = 'time hru'
#  ### NM added end
  
  # Write basic global attribute
  ncfile.history = 'Created ' + time.ctime(time.time())
  ncfile.source = os.path.dirname(os.path.abspath(__file__))+__file__[1:]
  
  ncfile.close()

def compAvgVal(nc_wgt,nc_in,varname):
  """Compute areal weighted avg value of <varname> in <nc_in> for each hru based on hru's weight in <nc_wgt>""" 
  wgt = wgtnc(nc_wgt) #instantaneous of wgtnc object
  hruIdName = wgt.getHruIdName()
  hruIDs    = wgt.getHruID() #get hruID list 

  dataVal  = getNetCDFData(nc_in,varname) # Get data value 
  FillVal  = getNetCDFAtt(nc_in,varname,'_FillValue') # Get data value 
  IdVal    = getNetCDFData(nc_in,'hru_id') 
  dim1size = dataVal.shape[0]
  #latVal   = getNetCDFData(nc_in,'lat') 
  #lonVal   = getNetCDFData(nc_in,'lon') 
  #if lonVal[1]>0:
  #  lonVal = lonVal-360

  #Initialize array
  wgtVal = np.zeros((dim1size,len(hruIDs)))
  #Loop through each hru polygon 
  for i in range(len(hruIDs)):
    print "Computing weighted values over hru%d" %hruIDs[i]
    # Get list of wgt, lat, and lon for corresponding hru
    (wgtval, overlapsId, lat, lon, overlaps)=wgt.getWgtHru(hruIDs[i])
    wgtArray = np.asarray(wgtval)
    
    #Count missing cells (No values in original input polygon) and valid cells
    numvoid = 0
    # Go through wgt list and replace value with zero where overlapping polygon has missing value - case1
    #                                              or where overlapping polygon is outside nc_wgt domain - case2
    for j in range(len(wgtArray)): # Go through each overlapping polygon

      if overlapsId[j] >= 0: # if lat list is not empty i.e., there is at least one overlapping polygon
        #ij=np.column_stack(np.where(IdVal==overlapsId[j]))    
        ij=np.where(IdVal==overlapsId[j])
        if len(ij[0]) == 0: # if nc_in netCDF does not cover nc_wgt domain - case2  old statement  if len(ij) == 0:
          numvoid = numvoid+1
          wgtArray[j]=0
        else:
          #a = dataVal[10,ij[0,0],ij[0,1]].tolist()
          a = dataVal[10,ij[0][0],ij[1][0]].tolist()

          if a == None: #if value of overlapping polygon is missing data -case1
            numvoid = numvoid+1
            wgtArray[j] = 0 
    
      #if lat[j]: # if lat list is not empty i.e., there is at least one overlapping polygon
      #  yj=np.where(latVal==lat[j])          
      #  xi=np.where(lonVal==lon[j])          
      #  if len(yj[0]) == 0 or len(xi[0]) == 0: # if nc_in netCDF does not cover nc_wgt domain - case2
      #    numvoid = numvoid+1
      #    wgtArray[j]=0
      #  else:
      #    a = dataVal[0,yj,xi].tolist()
      #    if a[0][0] == None: #if value of overlapping polygon is missing data -case1
      #      numvoid = numvoid+1
      #      wgtArray[j] = 0 
    
    #Adjust weight value if valid weight value (> 0) exist in list
    #if ( np.nansum(wgtArray) != 0.0 and overlaps != 0.0 ):

    Val2 = np.zeros((dim1size))
    if np.nansum(wgtArray) > 0.0:
    #  print (np.nansum(wgtArray) > 0.0)
    #  print (overlaps != 0.0)
      newWgtArray = [x/np.nansum(wgtArray) for x in wgtArray]
    #  print "%d invalid polygons out of %d polygons" %(numvoid,overlaps)
    #  print "sum of original weights is %f" %np.nansum(wgtArray)
    #  print "sum of new weights is %f" %np.nansum(newWgtArray)
    #  print newWgtArray

      #Initialize data storage for time series of weighted value for input polygon 
      c=0
      #Val1 = np.zeros((dim1size,1))
      for j in range(len(newWgtArray)): # Go through each overlapping polygon
        if overlapsId[j] >= 0: # if lat list is not empty i.e., there is at least one overlapping polygon

          #ij=np.column_stack(np.where(IdVal==overlapsId[j]))    
          ij=np.where(IdVal==overlapsId[j])    
          yj = ij[0]
          xi = ij[1]
          Val1=newWgtArray[j]*dataVal[:,yj,xi]
          #If selected hru id from runoff nc does not match with wieght nc, exit 
          if (IdVal[yj,xi] != overlapsId[j]):
             print "hru_id %d not exist, Sayonara" %IdVal[yj,xi]
             sys.exit()

          if Val1.size: #if Val1 is not empty i.e., exclude empty value in weight list
            c += 1
            if c == 1:
              oriVal = Val1
            else:
              oriVal=np.concatenate((oriVal, Val1), axis=1)

      Val2 = oriVal.sum(axis=1)
      del oriVal
      #  if lat[j]: # if lat list is not empty i.e., there is at least one overlapping polygon
      #    yj=np.where(latVal==lat[j])          
      #    xi=np.where(lonVal==lon[j])          
      #    Val1=newWgtArray[j]*dataVal[:,yj,xi]
      #    if Val1.size: #if Val1 is not empty i.e., exclude empty value in weight list
      #      c += 1
      #      if c == 1:
      #        oriVal = Val1
      #      else:
      #        oriVal=np.concatenate((oriVal, Val1), axis=1)
      #Val2=oriVal.sum(axis=1)

    #if there is no valid weight values (i.e. all zero)
    else:
      Val2 = np.ones((dim1size))*-9999.0

    #if i == 0:
    #  wgtVal = Val2
    #else:
      #wgtVal=np.concatenate((wgtVal, Val2), axis=1)
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
      nc_wgt = args[0]
      nc_in = args[1]  
      varname = args[2]  
      nc_out = args[3]

      wgtVal=compAvgVal(nc_wgt,nc_in,varname)
      writeNetCDFData(nc_out, wgtVal, varname)     
    else:
      usage()

