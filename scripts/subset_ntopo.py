#!/usr/bin/env python

import sys
import getopt
import xarray as xr
import numpy as np

# -----
# Hardcoded variables
# -----
# Fill values
FILL_VALUE = -9999.0

# Variables with sSeg dimension
seg_var_name = ['reachIndex',
                'reachID',
                'reachSlope',
                'reachLength',
                'basinArea',
                'upstreamArea',
                'totalArea',
                'downReachIndex',
                'downReachID',
                'reachLat1',
                'reachLat2',
                'reachLon1',
                'reachLon2',
                'reachStart',
                'reachCount',
                'upReachStart',
                'upReachCount',
                'upHruStart',
                'upHruCount']
# Variables with sHRU dimension
hru_var_name = ['hruIndex',
                'hru_id',
                'hru_lon',
                'hru_lat',
                'hru_elev',
                'hru_area',
                'hru_weight']
# Variables with sUps dimension
ups_var_name = ['upReachIndex',
                'upReachID']
# Variables with sAll dimension
allups_var_name = ['reachList',
                   'upReachTotalLength']


def read_ntopo(nc_in):
    ntopo_data = xr.open_dataset(nc_in)
    return ntopo_data


def subset_ntopo(seg_id, in_data):
    """Find upstream segs of <seg_id> in network data <in_data>
       extract inforation of upstream segs and return them in <out_data>
    """
    print('Outlet segment = %s' % seg_id)

    # Identify index of the desired stream segment from reachID vector
    reach_id_array = in_data['reachID'].values     # reach ID
    idx = np.where(reach_id_array == int(seg_id))
    if not idx:
        raise Exception('unable to find desired stream segment')

    # Split into separate datasets based on a dimension
    seg_data = in_data
    hru_data = in_data
    ups_data = in_data
    for varname in in_data.data_vars:
        if in_data[varname].dims[0] != 'sSeg':
            seg_data = seg_data.drop(varname)
        if in_data[varname].dims[0] != 'sHRU':
            hru_data = hru_data.drop(varname)
        if in_data[varname].dims[0] != 'sUps':
            ups_data = ups_data.drop(varname)

    idx_rch_start = in_data['reachStart'].values[idx]-1
    n_rch_count = in_data['reachCount'].values[idx]
    print 'Number of upstream segment from outlet segment (nRchCount) = %d' % n_rch_count

    # Read reach list of index from global segments (all the upstream reachs for each segment)
    idxs_uprch_list = in_data['reachList'].values[idx_rch_start[0]:idx_rch_start[0]+n_rch_count[0]]

    # Reach upstream segment and associated HRU infor from non-ragged vector
    seg_sub_data = seg_data.isel(sSeg = idxs_uprch_list-1)
    print 'upstream seg id list = ',seg_sub_data['reachID'].values

    # Redo index for 'reachIndx' as local segment list
    seg_sub_data['reachIndex'] = xr.DataArray(np.arange(1,n_rch_count+1), dims=['sSeg'])

    # Re-compute downstream segment index as local segment list
    down_idx_list = np.full(n_rch_count, FILL_VALUE, dtype=np.int32)
    for ii, down_id in enumerate(seg_sub_data['downReachID'].values):
        down_idx = np.where(seg_sub_data['reachID'].values == down_id)[0]
        if np.size(down_idx) == 1:
            down_idx_list[ii] = down_idx+1
    seg_sub_data['downReachIndex'] = xr.DataArray(down_idx_list, dims=['sSeg'])

    # Assign downstream segment ID = 0 at desired outlet segment
    for iSeg, id_uprch in enumerate(seg_sub_data['reachID'].values):
        if (id_uprch == int(seg_id)):
            seg_sub_data.data_vars['downReachID'][iSeg] = 0
            seg_sub_data.data_vars['downReachIndex'][iSeg] = -9999

    # Reach upstream segment and associated HRU infor from ragged vector
    uprch_len_list = []
    uprch_idx_list = []
    hru_data_list = []
    ups_data_list = []
    for idx_uprch in idxs_uprch_list:
        idx_uprch -= 1
        idx_rch_start_tmp = in_data['reachStart'].values[idx_uprch]-1
        n_rch_count_tmp = in_data['reachCount'].values[idx_uprch]

        # sAll dimension
        # upReachTotalLength
        uprch_len_list.append(in_data['upReachTotalLength'][idx_rch_start_tmp:idx_rch_start_tmp+n_rch_count_tmp])

        # reachList
        reach_list_tmp = in_data['reachList'].values[idx_rch_start_tmp:idx_rch_start_tmp+n_rch_count_tmp]
        up_idx_array = np.full(n_rch_count_tmp, FILL_VALUE, dtype=np.int32)
        for ii, up_id in enumerate(in_data['reachID'].values[reach_list_tmp-1]):
            up_idx = np.where(seg_sub_data['reachID'].values == up_id)[0]
            up_idx_array[ii] = up_idx+1
        uprch_idx_list.append(xr.DataArray(up_idx_array, dims=['sAll'], name='reachList'))

        # Recompute all the upstream segment indices as local segment list
        # sUps dimension
        idx_ups_start = in_data['upReachStart'].values[idx_uprch]-1
        n_ups_count = in_data['upReachCount'].values[idx_uprch]
        if (n_ups_count > 0):
             ups_data_list.append(ups_data.isel(sUps = range(idx_ups_start, idx_ups_start+n_ups_count)))

        # sHrus dimension
        idx_uphru_start = in_data['upHruStart'].values[idx_uprch]-1
        n_uphru_count = in_data['upHruCount'].values[idx_uprch]
        if (n_uphru_count > 0):
            hru_data_list.append(hru_data.isel(sHRU = range(idx_uphru_start, idx_uphru_start+n_uphru_count)))

    cum_n = 1
    for ii, n_uprch in enumerate(seg_sub_data['reachCount'].values):
        seg_sub_data.data_vars['reachStart'][ii] = cum_n
        cum_n += n_uprch

    cum_n = 1
    for ii, n_uphru in enumerate(seg_sub_data['upHruCount'].values):
        seg_sub_data.data_vars['upHruStart'][ii] = cum_n
        cum_n += n_uphru

    cum_n = 1
    for ii, n_iuprch in enumerate(seg_sub_data['upReachCount'].values):
        if n_iuprch == 0:
            seg_sub_data.data_vars['upReachStart'][ii] = FILL_VALUE
        else:
            seg_sub_data.data_vars['upReachStart'][ii] = cum_n
        cum_n += n_iuprch

    up_len_data = xr.concat(uprch_len_list, dim='sAll')
    up_idx_data = xr.concat(uprch_idx_list, dim='sAll')
    hru_sub_data = xr.concat(hru_data_list, dim='sHRU')
    ups_sub_data = xr.concat(ups_data_list, dim='sUps')

    # Re-compute immediate upstream segment index as local segment list
    imm_up_idx_list = np.full(ups_sub_data['upReachID'].size, FILL_VALUE, dtype=np.int32)
    for ii, imm_up_id in enumerate(ups_sub_data['upReachID'].values):
        imm_up_idx = np.where(seg_sub_data['reachID'].values == imm_up_id)[0]
        if np.size(imm_up_idx) == 1:
            imm_up_idx_list[ii] = imm_up_idx+1
    ups_sub_data['upReachIndex'] = xr.DataArray(imm_up_idx_list, dims=['sUps'])

    out_data = xr.merge([seg_sub_data, ups_sub_data, hru_sub_data, up_len_data, up_idx_data])

    return out_data


# Main
use = '''
Usage: %s -[h] <notpo_in_netCDF>
               <seg_id>
               <ntopo_out_netCDF>
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

    if len(args) == 3:
        nc_in = args[0]
        seg_id = args[1]
        nc_out = args[2]

        ntopo_data = read_ntopo(nc_in)
        sub_ntopo_data = subset_ntopo(seg_id, ntopo_data)
        print(sub_ntopo_data)
        sub_ntopo_data.to_netcdf(nc_out, format='NETCDF4')

    else:
        usage()
