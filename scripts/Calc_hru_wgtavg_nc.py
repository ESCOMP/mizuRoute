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
NPROC=2
# Name of variables in runoff netCDF
GRID_ID_NAME='hru_id'
TIME_VAR_NAME='time'
# Name of variables in weight netCDF
HRU_ID_NAME='hruid'
WGTNM='weight'
OVRPLYNM='overlapPolyId'
OVRNM='overlaps'

FILL_VALUE=-9999.0
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
        wgtAll        = get_netCDF_data(self.ncName, WGTNM) 
        overlapsIdAll = get_netCDF_data(self.ncName, OVRPLYNM) 

        self.hruList = self.getHruID()
        idx=self.hruList.index(hru)               # indix of corresponding hru
        self.wgt = np.asarray(list(wgtAll[idx]))       # list of areal weight of overlapping poly for hru
        self.overlap_ids = np.asarray(list(overlapsIdAll[idx])) # list of IDs of overlapping poly for hru
        self.overlap_ids=self.overlap_ids[~np.isnan(self.overlap_ids)]
        self.wgt=self.wgt[~np.isnan(self.wgt)]

        return (self.wgt, self.overlap_ids)

    def getHruID(self):
        """ get hru ID list of basin"""
        self.hruid = list(get_netCDF_data(self.ncName, HRU_ID_NAME))
        return self.hruid

############################################
#            Modules                       #
############################################
def get_netCDF_data(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
    return data

def get_netCDF_attr(fn, varname,attName):
    """Read attribute of <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    var = f.variables[varname]
    attData = getattr(var,attName)
    return attData

def comp_agv_val(nc_wgt, nc_in, varname, chunk):
    """Compute areal weighted avg value of <varname> in <nc_in> for each hru based on hru's weight in <nc_wgt>""" 
    wgt       = wgtnc(nc_wgt)
    hruIDs    = wgt.getHruID()
    hrudata   = split(hruIDs,NPROC)
    hrulist   = hrudata[chunk]

    dataVal  = get_netCDF_data(nc_in,varname)
    FillVal  = get_netCDF_attr(nc_in,varname,'_FillValue')
    ghruid   = get_netCDF_data(nc_in,GRID_ID_NAME)
    dim1size = dataVal.shape[0]
    #Initialize wgt_ave_val[ntime,nhru] 
    wgt_ave_val = np.full((dim1size,len(hrulist)),FILL_VALUE)
    for i in range(len(hrulist)):
        print hrulist[i] # Get list of wgt for corresponding hru
        (wgt_array, overlap_ids)=wgt.getWgtHru(hrulist[i])
        # Go through wgt list and replace value with zero for following cases
        # where overlapping polygon has missing value - case1
        # where overlapping polygon is outside nc_wgt domain - case2
        a = np.zeros((dim1size,len(wgt_array)))
        for j,overlap_id in enumerate(overlap_ids):
            # find index of grid cell that match up with hru id of overlap_ids
            row,col=np.where(ghruid==overlap_id)
            # if nc_in netCDF does not cover nc_wgt domain - case2
            if not (np.size(row) and np.size(col)):
                wgt_array[j]=0.0
            else:
                a[:,j] = np.squeeze(dataVal[:,row,col])
                #if value of overlapping polygon is missing data -case1
                if any(a[:,j]==FillVal):
                    wgt_array[j] = 0.0
        sum_wgt=np.sum(wgt_array)
        if sum_wgt > 0.0:
            # Adjust weight value if valid weight value (> 0) exist in list
            wgt_array = wgt_array/sum_wgt
            wgt_array = np.where(np.isnan(wgt_array),0,wgt_array)
            wgt_ave_val[:,i] = np.dot(a,wgt_array.T).T

    return wgt_ave_val

def split(hrulist,n):
    subsize=len(hrulist)/n
    idx1=0
    hrudata=[]
    for _ in range(n-1):
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
    for (opt,val) in opts:
        if opt == '-h':
            usage()
        elif opt == '-v':
            verbose = True
        else:
            raise OptionError, opt
            usage()

    if len(args) == 4:
        nc_wgt = args[0]
        nc_in = args[1]  
        varname = args[2]  
        nc_out = args[3]

        func=partial(comp_agv_val, nc_wgt, nc_in, varname)

        pool = mp.Pool(processes=NPROC,maxtasksperchild=1)
        results = pool.map(func,range(NPROC))
        result = np.column_stack(results)

        with xr.open_dataset(nc_in) as ds:
            timedata = ds[TIME_VAR_NAME]
        with xr.open_dataset(nc_wgt) as ds:
            hrudata = ds[HRU_ID_NAME]

        encoding={varname: {'dtype':'float32', '_FillValue':FILL_VALUE}}
        foo = xr.DataArray(result, coords=[('time',timedata), ('hru',hrudata)], name=varname)
        foo.encoding = encoding[varname]
        foo.to_netcdf(nc_out)
    else:
        usage()

