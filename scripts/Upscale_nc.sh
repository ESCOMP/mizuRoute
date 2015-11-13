#!/bin/bash

# Make sure dimension in weight nc is hru_id2, otherwise edit it

if [ $# -ne 6 ]; then
  echo "This script compute weighted averaged time series for hru & add attributes in netCDF"
  echo "Not enough argument supplied" 
  echo "usage: ./Upscale_nc.sh <nc_wgt> <hrudim> <nc_in> <var> <nc-out> "
  echo "<nc_wgt>: name of weight netCDF"
  echo "<hrudim>: name of hru dimension in weight netCDF"
  echo "<nc_in>:  Name of input netCDF "
  echo "<var>:    Variable name in input netCDF e.g., RUNOFF"
  echo "<unit>:   Unit of variablee.g., mm/day"
  echo "<nc_out>: Name of output netCDF containing weighted avg values" 
  exit
fi

nc_wgt=$1
hrudim=$2
nc_in=$3
var=$4
unit=$5
nc_out=$6

#temporal netCDF
nc_time=${var}_time_temp.nc
nc_hru=${var}_hru_temp.nc

ncks -O -v time ${nc_in} ${nc_time}
ncks -O -v ${hrudim} ${nc_wgt} ${nc_hru}

echo 'Your arguments are like these'
echo ${nc_wgt}
echo ${nc_in}
echo ${var}
echo ${nc_out}

#Run python script to compute upscaled values
./Calc_hru_wgtavg_nc.py ${nc_wgt} ${nc_in} ${var} ${nc_out}

#Some post process
ncrename -O -d dim1,time ${nc_out} 
ncrename -O -d dim2,${hrudim} ${nc_out}
ncks -A ${nc_time} ${nc_out}
#ncatted -O -a units,${VAR},c,c,"mm/s" ${nc_out}  
ncatted -O -a units,${var},c,c,${unit} ${nc_out}  
ncatted -O -a _FillValue,${var},c,f,-9999.0 ${nc_out}
ncks -O --mk_rec_dmn time ${nc_out} ${nc_out}
ncks -A ${nc_hru} ${nc_out} 

rm -f ${nc_time} ${nc_hru} 

echo "Finished" 
