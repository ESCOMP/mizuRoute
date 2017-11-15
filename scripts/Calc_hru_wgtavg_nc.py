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
import pandas as pd
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


def get_wgt(nc_name, hru):
    """For given hru id, get weights of
       overlapping polygons and their IDs
    """
    wgt_all = get_netCDF_data(nc_name, WGTNM)
    overlap_id_all = get_netCDF_data(nc_name, OVRPLYNM)

    hru_list = get_hru_id(nc_name)
    idx = hru_list.index(hru)
    # list of areal weight of overlapping polygons
    wgt = np.asarray(list(wgt_all[idx]))
    # list of IDs of overlapping polygons
    overlap_id = np.asarray(list(overlap_id_all[idx]))
    overlap_id = overlap_id[~np.isnan(overlap_id)]
    wgt = wgt[~np.isnan(wgt)]
    return (wgt, overlap_id)


def get_hru_id(nc_name):
    """ get hru ID list of basin"""
    hru_id = get_netCDF_data(nc_name, HRU_ID_NAME).tolist()
    return hru_id


def get_netCDF_data(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn, 'r')
    var_data = f.variables[varname][:]
    f.close()
    return var_data


def get_netCDF_attr(fn, varname, att_name):
    """Read attribute of <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn, 'r')
    var = f.variables[varname]
    att_data = getattr(var, att_name)
    return att_data


def comp_agv_val(nc_wgt, nc_in, varname, chunk):
    """Compute areal weighted avg value of <varname> in <nc_in>
       for each of <chunk> hrus based on hru's weight in <nc_wgt>
    """
    hru_ids = get_hru_id(nc_wgt)
    hru_data = split(hru_ids, NPROC)
    hru_list = hru_data[chunk]

    var_data = get_netCDF_data(nc_in, varname)
    fill_val = get_netCDF_attr(nc_in, varname, '_FillValue')
    ghruid = get_netCDF_data(nc_in, GRID_ID_NAME)
    dim1size = var_data.shape[0]
    # Initialize wgt_ave_val[ntime,nhru]
    wgt_ave_val = np.full((dim1size, len(hru_list)), FILL_VALUE)
    for i in range(len(hru_list)):
        print(hru_list[i])  # Get list of wgt for corresponding hru
        wgt_array, overlap_ids = get_wgt(nc_wgt, hru_list[i])
        # Go through wgt list and replace value with zero for following cases
        # where overlapping polygon has missing value - case1
        # where overlapping polygon is outside nc_wgt domain - case2
        a = np.zeros((dim1size, len(wgt_array)))
        for j, overlap_id in enumerate(overlap_ids):
            # find index of grid cell that match up with hru id of overlap_ids
            row, col = np.where(ghruid == overlap_id)
            # if nc_in netCDF does not cover nc_wgt domain - case2
            if not (np.size(row) and np.size(col)):
                wgt_array[j] = 0.0
            else:
                a[:, j] = np.squeeze(var_data[:, row, col])
                # if value of overlapping polygon is missing data -case1
                if fill_val in a[:, j]:
                    wgt_array[j] = 0.0
        sum_wgt = wgt_array.sum()
        if sum_wgt > 0.0:
            # Adjust weight value if valid weight value (> 0) exist in list
            wgt_array = wgt_array/sum_wgt
            wgt_array = np.where(pd.isnull(wgt_array), 0, wgt_array)
            wgt_ave_val[:, i] = np.dot(a, wgt_array.T).T

    return wgt_ave_val


def split(hru_list, n):
    subsize = len(hru_list)/n
    idx1 = 0
    hru_data = []
    for _ in range(n-1):
        idx2 = idx1 + subsize
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
        opts, args = getopt.getopt(sys.argv[1:], 'h')
    except getopt.error:
        usage()

    verbose = False
    for opt, val in opts:
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
