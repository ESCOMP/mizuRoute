function [var,attvar1,attvar2,attvar3,attvar4] = netcdf2mat(filename, varname, varargin)
% Project : General
%
% Purpose    : Read variable in netcdf
% Syntax     : netcdf2mat(filename, varname)
% Arguments  : filename    = string - input netCDF file name
%              narname     = string - variable name which to be read 
% Usage   : [var,time, long, lat]=...
%                read_wrf(infile, ncvarname, 'time',1, 'latlon',1)
%          
%           infile        = input file name including full path (string).
%           ncvarname     = variable names in NetCDF (string)
%           time          = option to get time (1-> yes, 0 no)
%           latlon        = option to get lat/long grid (1-> yes, 0 no) 
%           the following example output time and latlon arrays
%           e.g., [var,t,x,y]=...
%                 read_wrf('/glade/home/mizukami/RAINNC_wrfout_d01_2008-12-31_00:00:00.nc', 'rain_rate', 'time',1, 'latlon',1)

% Description: Return axid of map so that map can be modified
%             e.g., netcdf2mat('/glade/proj2/ral/RHAP/asd000/files02km/wrfout_coord_2km.nc', 'HGT')
%
% Last Modified : 1/5/2012  Naoki Mizukami

%if no argument, run ncdump to get variable name stored in netCDF
 if nargin == 1
    unix(['ncdump -h ' filename])
    return
 elseif nargin < 1
     error('Not enough argument');
 end

% Validate arguments
%see http://www.mathworks.com/help/techdoc/matlab_prog/bresuxt-11.html
p = inputParser;

% Define required input arguments - infile, variablename
p.addRequired('filename',@(x)ischar(x));
p.addRequired('varname',@(x)ischar(x));

%Define Parameter/value inputs (optional arguments)
p.addParamValue('attname1', 'na', @(x)ischar(x));
p.addParamValue('attname2', 'na', @(x)iscahr(x));
p.addParamValue('attname3', 'na', @(x)iscahr(x));
p.addParamValue('attname4', 'na', @(x)iscahr(x));

% Parse and validate all input arguments.
p.parse(filename,varname,varargin{:});

nc = netcdf.open(filename, 'NC_NOWRITE');

varid = netcdf.inqVarID(nc, varname);
var  = netcdf.getVar(nc,varid,'double');
%if 2D array (e.g. latitude, longitude, elevation grid)
if ndims(var)==2
    var  = permute( var, [2,1]);
% 3D array (e.g., time series of 2D grids like precipitation-x,y)
elseif ndims(var)==3
    var = permute( var, [2,1,3]);
% 4D array (e.g., time seris of 3D grids like soil moisture- x,y,depth)
elseif ndims(var)==4
    var = permute( var, [2,1,3,4]);
end
%attribute
 if ~strcmp(p.Results.attname1,'na')
     attvar1=netcdf.getAtt(nc,varid,p.Results.attname1);
 end
 if ~strcmp(p.Results.attname2,'na')
     attvar2=netcdf.getAtt(nc,varid,p.Results.attname2);
 end
  if ~strcmp(p.Results.attname3,'na')
     attvar3=netcdf.getAtt(nc,varid,p.Results.attname3);
 end
 if ~strcmp(p.Results.attname4,'na')
     attvar4=netcdf.getAtt(nc,varid,p.Results.attname4);
 end

netcdf.close(nc);

return
 
