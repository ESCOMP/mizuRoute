#!/usr/bin/env python

# this is NOT meant to be anything more than a start for people.
# it is NOT a general use utility.
# no one involved in creating or distributing this script has any
# responsibility for its function or output.

import sys
import multiprocessing as mp
from functools import partial
import getopt
import xarray as xr
import numpy as np
import netCDF4 as nc4

# -----
# Hardcoded variables
# -----
# number of processors (= number of chunks of hru to be processed)
NPROC = 2
# Name of variables in runoff netCDF
GRID_ID_NAME = 'hru_id'
TIME_VAR_NAME = 'time'
# Name of variables in weight netCDF
HRU_ID_NAME = 'hruid'
WGTNM = 'weight'
OVRPLYNM = 'overlapPolyId'
OVRNM = 'overlaps'

FILL_VALUE = -9999.0


class wgtnc:
    """ object of basin netCDF including
        areal weights and hrus of overlapping polygons
    """
    def __init__(self, ncName):
        """Initialization """
        self.ncName = ncName

    def getWgtHru(self, hru):
        """For given hru id, get weights of
           overlapping polygons and their IDs
        """
        wgtAll = get_netCDF_data(self.ncName, WGTNM)
        overlapsIdAll = get_netCDF_data(self.ncName, OVRPLYNM)

        self.hruList = self.getHruID()
        idx = self.hruList.index(hru)  # indix of current hru
        # list of areal weight of overlapping polygons
        self.wgt = np.asarray(list(wgtAll[idx]))
        # list of IDs of overlapping polygons
        self.overlap_ids = np.asarray(list(overlapsIdAll[idx]))
        self.overlap_ids = self.overlap_ids[~np.isnan(self.overlap_ids)]
        self.wgt = self.wgt[~np.isnan(self.wgt)]

        return (self.wgt, self.overlap_ids)

    def getHruID(self):
        """ get hru ID list of basin"""
        self.hruid = list(get_netCDF_data(self.ncName, HRU_ID_NAME))
        return self.hruid


def get_netCDF_data(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    fid = nc4.Dataset(fn, 'r')
    data = fid.variables[varname][:]
    fid.close()
    return data


def get_netCDF_attr(fn, varname, attName):
    """Read attribute of <varname> variables from NetCDF <fn> """
    fid = nc4.Dataset(fn, 'r')
    var = fid.variables[varname]
    att_data = getattr(var, attName)
    return att_data


def comp_agv_val(nc_wgt, nc_in, varname, chunk):
    """Compute areal weighted avg value of <varname> in <nc_in>
       for each of <chunk> hrus based on hru's weight in <nc_wgt>
    """
    wgt = wgtnc(nc_wgt)
    hru_ids = wgt.getHruID()
    hru_data = split(hru_ids, NPROC)
    hru_list = hru_data[chunk]

    data_val = get_netCDF_data(nc_in, varname)
    fill_val = get_netCDF_attr(nc_in, varname, '_FillValue')
    ghruid = get_netCDF_data(nc_in, GRID_ID_NAME)
    dim1size = data_val.shape[0]
    # Initialize wgt_ave_val[ntime,nhru]
    wgt_ave_val = np.full((dim1size, len(hru_list)), FILL_VALUE)
    for i in range(len(hru_list)):
        print hru_list[i]  # Get list of wgt for corresponding hru
        (wgt_array, overlap_ids) = wgt.getWgtHru(hru_list[i])
        # Go through wgt list and replace value with zero for following cases
        # where overlapping polygon has missing value - case1
        # where overlapping polygon is outside nc_wgt domain - case2
        sub_data = np.zeros((dim1size, len(wgt_array)))
        for j, overlap_id in enumerate(overlap_ids):
            # find index of grid cell that match up with hru id of overlap_ids
            row, col = np.where(ghruid == overlap_id)
            # if nc_in netCDF does not cover nc_wgt domain - case2
            if not (np.size(row) and np.size(col)):
                wgt_array[j] = 0.0
            else:
                sub_data[:, j] = np.squeeze(data_val[:, row, col])
                # if value of overlapping polygon is missing data -case1
                if any(sub_data[:, j] == fill_val):
                    wgt_array[j] = 0.0
        sum_wgt = np.sum(wgt_array)
        if sum_wgt > 0.0:
            # Adjust weight value if valid weight value (> 0) exist in list
            wgt_array = wgt_array/sum_wgt
            wgt_array = np.where(np.isnan(wgt_array), 0, wgt_array)
            wgt_ave_val[:, i] = np.dot(sub_data, wgt_array.T).T

    return wgt_ave_val


def split(hru_list, n):
    subsize = len(hru_list)/n
    idx1 = 0
    hru_data = []
    for _ in range(n-1):
        idx2 = idx1+subsize
        hru_data.append(hru_list[idx1:idx2])
        idx1 = idx2
    hru_data.append(hru_list[idx1:])
    return hru_data


# ----
# Main
# ----
use = '''
Usage: %s -[h] <weight_netCDF>
               <input_netCDF>
               <variable_name_in_input_netCDF>
               <output_netCDF>
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
    for (opt, val) in opts:
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

        func = partial(comp_agv_val, nc_wgt, nc_in, varname)

        pool = mp.Pool(processes=NPROC, maxtasksperchild=1)
        results = pool.map(func, range(NPROC))
        result = np.column_stack(results)

        with xr.open_dataset(nc_in) as ds:
            timedata = ds[TIME_VAR_NAME]
        with xr.open_dataset(nc_wgt) as ds:
            hru_data = ds[HRU_ID_NAME]

        encoding = {varname: {'dtype': 'float32', '_FillValue': FILL_VALUE},
                    TIME_VAR_NAME: {'dtype': 'int32'}}
        foo = xr.DataArray(result,
                           coords=[(TIME_VAR_NAME, timedata), ('hru', hru_data)],
                           name=varname)
        foo.encoding = encoding[varname]
        foo[TIME_VAR_NAME].encoding = encoding[TIME_VAR_NAME]
        foo.to_netcdf(nc_out)
    else:
        usage()
