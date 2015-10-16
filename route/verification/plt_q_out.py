#!/usr/bin/env python

#plot time series of mizuRoute output
#it is NOT a general use utility.

import sys
import os
import time 
import getopt
import matplotlib.pyplot as plt 
import matplotlib as mpl 
import numpy as np
from scipy.stats import nanmean
import netCDF4 as nc4

############################################
# Hardcoded setting
############################################
SEGID=14000645
VARNM=['KWTroutedRunoff','IRFroutedRunoff'] 
############################################
#            Modules                       #
############################################
def getNetCDFData(fn, varname):
    """Read <varname> variables from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
    return data
############################################
#                Main                      #
############################################
use = '''
Usage: %s -[h] -i<name_netCDF> 
        -h  help
        -i  netCDF output from mizuRoute
'''

if __name__ == '__main__':

  def usage():
      sys.stderr.write(use % sys.argv[0])
      sys.exit(1)
  try:
      (opts, args) = getopt.getopt(sys.argv[1:], 'hi:')
  except getopt.error:
      usage()
  
  innc=None
  for (opt,val) in opts:  
    if opt == '-h':
      usage()
    elif opt == '-i':
      innc = val
    else:
      usage()
  if not innc:
    usage()
  
  # Get index of desired segID
  reachID = getNetCDFData(innc,'reachID') # Get data value from hru netCDF 
  idx=np.where(reachID==SEGID)

  #read netCDF file 
  ro=dict()
  for v in VARNM:
    ro[v] = getNetCDFData(innc,v) # Get data value from hru netCDF 
  
  # plot time series
  mpl.rcParams['axes.color_cycle']=['r','b','g','c','m','k']
  for key,val in ro.items():
    # val.shape[0]: time dimension
    # val.shape[1]: seg dimension
    plt.plot(val[:,idx[0]],
                linewidth=1.0,
                linestyle="-",
                label=key)
    plt.xlim(0,366)
    plt.xlabel('Day from October 1')
    plt.ylabel('q (cms)')
    plt.title('segID='+str(SEGID))

    # Now add the legend with some customizations.
    plt.legend(loc='upper right', shadow=True)

  plt.show()
