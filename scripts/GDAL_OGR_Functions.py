#!/usr/bin/env python
#-------------------------------------------------------------------------------
# This is collection of GIS manipulation functions 
# Use the OGR library
#-------------------------------------------------------------------------------

import os, sys, ogr, osr, gdal
import shapely
from shapely.wkb import loads
from shapely import speedups
from shapely.ops import cascaded_union
from gdalconst import *
import time
from netCDF4 import Dataset
import numpy as np
import ogr2ogr

# Setup directories to import ogr2ogr.py
pwd = os.path.dirname(__file__)                                                 # Store current working directory
print 'path to files: %s' %pwd
sys.path.append(pwd)                                                            # Append current directory to the python path
if shapely.speedups.available:
    shapely.speedups.enable()

driver = ogr.GetDriverByName('ESRI Shapefile')

def featureIsValid(shgeom):
    """ Return True if this feature passes validity tests. """
    if not shgeom.is_valid:
        return "Geometry is not valid."
    if shgeom.is_empty:
        return "Geometry is empty."
    return None

def polygonize_raster(inraster, outputFile):
    '''Function takes a directory and a raster (as defined by GDAL) and creates
    a polygon shapefile with one feature per raster cell.  Currently only EPSG:4326
    is a supported input raster coordinate system.'''

    gdal.AllRegister()

    # Opening the file with GDAL, with read only acces
    dataset = gdal.Open(inraster, GA_ReadOnly)

    # Getting raster dataset information
    #print 'Driver: %s: %s' %(dataset.GetDriver().ShortName,dataset.GetDriver().LongName)
    print 'Input Raster Size: %s x %s x %s' %(dataset.RasterXSize,dataset.RasterYSize,dataset.RasterCount)
    #print 'Projection: %s' %dataset.GetProjection()
    geotransform = dataset.GetGeoTransform()
    x0 = geotransform[0] # top left x
    y0 = geotransform[3] # top left y
    pwidth = geotransform[1] # pixel width
    pheight = geotransform[5] # pixel height is negative because it's measured from top
    if not geotransform is None:
        print 'Origin (x,y): %s,%s' %(x0,y0)
        print 'Pixel Size (x,y): %s,%s' %(pwidth,pheight)

     ########################### inserted on 11/06/14
    # Get raster values as array (must have Numpy 1.8.1 or greater)
    band = dataset.GetRasterBand(1)
    ndv = float(band.GetNoDataValue())
    DataType = band.DataType
    print 'NoData Value: %s' %ndv
    print 'Raster DataType: %s' %gdal.GetDataTypeName(DataType)
    #data = band.ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize)
    data = band.ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize).astype(np.float)
     ###########################

    # Now convert it to a shapefile with OGR
    datasource = driver.CreateDataSource(outputFile)
    #driver = ogr.GetDriverByName( 'Memory' )                                   # For in-memory vector file
    if datasource is None:
        print "Could not create "+outputFile+"\n"
        return False

    # Create the SpatialReference
    coordinateSystem = osr.SpatialReference()
    coordinateSystem.ImportFromEPSG(4326)                                           # WGS 84
    #coordinateSystem.ImportFromWkt(dataset.GetProjection())                         # Use projection from input raster

    # Create a new layer on the data source with the indicated name, coordinate system, geometry type
    layer = datasource.CreateLayer(outputFile, coordinateSystem, geom_type=ogr.wkbPolygon)
    if layer is None:
        print "Layer creation failed.\n"

    # Create a new field on a layer. Add one attribute
    layer.CreateField(ogr.FieldDefn('hru_id', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('lon_cen', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('lat_cen', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('CELLVALUE', ogr.OFTReal))

    # Fetch the schema information for this layer
    LayerDef = layer.GetLayerDefn()

    # For each cell, create a rectangle
    i = 1
    counter = 1
    for x in range(dataset.RasterXSize):
        j = 1
        for y in range(dataset.RasterYSize):
            # Calculating the polygon's coordinates that frame the raster image
            x00 = x0 + (pwidth * x)
            y00 = y0 - (abs(pheight) * y)
            x1 = x00 + pwidth
            y1 = y00
            x2 = x00
            y2 = y00 - abs(pheight)
            x3 = x1
            y3 = y2

            #create polygon object:
            myRing = ogr.Geometry(type=ogr.wkbLinearRing)
            myRing.AddPoint(x00, y00)
            myRing.AddPoint(x1, y1)
            myRing.AddPoint(x3, y3)
            myRing.AddPoint(x2, y2)
            myRing.AddPoint(x00, y00)                                           #close ring
            geometry = ogr.Geometry(type=ogr.wkbPolygon)
            geometry.AddGeometry(myRing)

            # create point geometry for coordinate system tranformation
            pt = geometry.Centroid()

            # Create a new feature (attribute and geometry)
            feature = ogr.Feature(LayerDef)
            feature.SetField('hru_id', counter)
            feature.SetField('lon_cen', pt.GetX())
            feature.SetField('lat_cen', pt.GetY())
            feature.SetField('CELLVALUE', data[y,x])                       # Add in the raster value

            # Make a geometry from Shapely object
            feature.SetGeometry(geometry)
            layer.CreateFeature(feature)

            #flush memory
            feature = geometry = myRing = None  # destroy these
            
            counter += 1
            j += 1
        i += 1

    # Save and close everything
    datasource = layer = None
    return True

def polyXY(inPolyShp,hrufieldname):
  '''
  Input: inPolySh: shapefile name (string)
       : polyID: Shapefile ID
  Output: XY {polyID: [lon lat] dictionary  
  '''
  # Open the basin shapefile file with GDAL, with read only access
  poly = driver.Open(inPolyShp,0)
  if poly is None:
      print "Open failed.\n"
  # Open a layer    
  layer = poly.GetLayer()

  # Getting features
  XY={}
  for feature in layer:
    FeatureId = feature.GetField(hrufieldname)
    # Getting Geometry object (point, polygon, etc.) 
    geometry = feature.GetGeometryRef()
    shgeom = shapely.wkb.loads(geometry.ExportToWkb())
    XY[FeatureId] = (shgeom.centroid.x,shgeom.centroid.y)

  poly.Destroy

  return XY 

def project_to_input(infeatures, outputFile2, gridpolys):
    shp = driver.Open(infeatures, 0)                                                # Opening the file with GDAL, with read only acces
    lyr = shp.GetLayer()
    spatialref = lyr.GetSpatialRef().ExportToWkt()
    if shp is None:
        print "Open failed.\n"
    ogr2ogr.main(["","-f", "ESRI Shapefile", "-t_srs", spatialref, outputFile2, gridpolys])   #note: main is expecting sys.argv, where the first argument is the script name, so the argument indices in the array need to be offset by 1
    shp = lyr = spatialref = None
    return outputFile2

def spatial_weights(hruShp,  hrufieldname, overlapPolyShp, hrufieldname2):
    '''This function takes an input hru shapefile and an input another polygon
    shapefile called "overlapping polygon" (both must be in the same coordinate reference sysetm) and performs
    spatial analysis to compute the individual weight of each overlapping polygon.
    This function uses only OGR routines
    Output:
    weight:        dictionary: {key:value}={hruid: [overlapPolyId, weight]}
    '''
    # time all the processes
    tic = time.time()

    # Initiate dictionaries
    weights = {}

    # Open the input shapefile file - shp1 and layer1
    shp1 = driver.Open(hruShp)
    if shp1 is None:
        print "Open failed.\n"
    layer1 = shp1.GetLayer()
    extent1 = layer1.GetExtent()

    # Open the overlaying polygon shapefile - shp2 and layer2
    shp2 = driver.Open(overlapPolyShp, 0)
    if shp2 is None:
        print "Open failed.\n"
    layer2 = shp2.GetLayer()

    # Create a Polygon from the extent tuple
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(extent1[0],extent1[2])
    ring.AddPoint(extent1[1],extent1[2])
    ring.AddPoint(extent1[1],extent1[3])
    ring.AddPoint(extent1[0],extent1[3])
    ring.AddPoint(extent1[0],extent1[2])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    layer2.SetSpatialFilter(poly)
    poly = None

    # Load entire overlaying polygon into shapely [id, geometry]
    overlayPolys = [(feat.GetField(hrufieldname2), loads(feat.GetGeometryRef().ExportToWkb())) for feat in layer2]

    # How many polygon in input shapefile?
    trials = layer1.GetFeatureCount()

    # Get information about field type
    layerDefinition = layer1.GetLayerDefn()
    for i in range(layerDefinition.GetFieldCount()):
        if layerDefinition.GetFieldDefn(i).GetName() == hrufieldname:
            field_defn = layerDefinition.GetFieldDefn(i)

    if field_defn.GetType() == ogr.OFTInteger:
        fieldtype = 'integer'                                                   #print "%d" % feat.GetFieldAsInteger(i)
    elif field_defn.GetType() == ogr.OFTReal:
        fieldtype = 'real'                                                      #print "%.3f" % feat.GetFieldAsDouble(i)
    elif field_defn.GetType() == ogr.OFTString:
        fieldtype = 'string'                                                    #print "%s" % feat.GetFieldAsString(i)
    else:
        fieldtype = 'string'                                                    #print "%s" % feat.GetFieldAsString(i)
    print "Field Type: %s (%s)" %(field_defn.GetType(), fieldtype)

    # Get list of unique fieldnames in input field
    hruidList = []
    feature = layer1.GetNextFeature()
    while feature is not None:
        hruidList.append(feature.GetField(hrufieldname))
        feature.Destroy()
        feature = layer1.GetNextFeature()
    uniqueHruIdList = list(set(hruidList))      # Creating a set eliminates duplicates
    layer1.ResetReading()
    print '    Finished gathering unique fieldnames...'

    # loop through the unique values of input features
    invalidcount = 0
    invalidlist = []
    for hruid in uniqueHruIdList:
        layer1.SetAttributeFilter("%s = '%s'" %(hrufieldname, hruid))
        basinpolys = [loads(feat.GetGeometryRef().ExportToWkb()) for feat in layer1]
        layer1.ResetReading()

        # Handle situations where only one goemetry exists
        if len(basinpolys) > 1:
            shgeom = shapely.ops.cascaded_union(basinpolys)
        else:
            shgeom = basinpolys[0]

        polygon_area = shgeom.area
        #basinareas = [basin.area for basin in basinpolys]
        #print ' Number of basins in cascaded union for gage %s: %s' %(value, len(basinpolys))
        #print ' Min: %s, Max: %s, Sum: %s, Discrepancy: %s' %(min(basinareas), max(basinareas), sum(basinareas), (sum(basinareas)-polygon_area))

        try:
            # Attempt to find all overlapping polygons [id, intersected area]
            weights[hruid] = [(poly[0],(poly[1].intersection(shgeom).area / polygon_area)) for poly in overlayPolys if shgeom.intersects(poly[1])==True]

        except:
            invalidcount += 1
            invalidlist.append([index, polygon_id])
            pass

    print "Invalid Count: %s" %invalidcount
    print "%s trials completed in %s seconds" %(trials, time.time()-tic)
    return weights, fieldtype

#def create_ncfile(lonlatdict, weights, outputfile, fieldtype): 
def create_ncfile(in_LatLon, in_weight, outputfile, hrufieldname, fieldtype): 
    '''This function creates an netcdf file with two dimensions (<hrufieldname> and 'overlapPoly')
    Input variable
    -in_LatLon         : dictionary: {key:value}= {overlayPolyId: [lon,lat]}
    -in_weights        : dictionary: {key:value}= {hruid: [overlapPolyId,weight]}
    -outputfile :
    -hrufieldname : hru dimension name
    -fieldtype  : type of hruid (str or int)

    NetCDF include
    -weight(hruid,overlapPoly)       : Weight value of each overlappin polygons for each hruid
    -latitude(hruid,overlapPoly)     : Centroid lat of each overlappin polygons for each hruid
    -longitude(hruid,overlapPoly)    : Centroid lon of each overlappin polygons for each hruid
    -overlapPolyId(hruid,overlapPoly): Overlapping polygon ID 
    -overlaps(hruid)                 : number of overlapping polygons

    which describe each basin and the grid cells that intersect it.  This leaves
    'whitespace' or blank index values and takes up unecessary disk space.'''

    # Create netcdf file for this simulation
    rootgrp = Dataset(outputfile, 'w', format='NETCDF4')

    # Get dimension size for overlapPoly dimension = length of largest grid to basin overlap
    overlapPolyNum = max([len(x) for x in in_weight.values()])

    # Create dimensions
    hruid_dim = rootgrp.createDimension(hrufieldname, len(in_weight))
    overlapPoly_dim = rootgrp.createDimension('overlapPoly', overlapPolyNum)

    # Create Coordinate variable
    # 1st dimension- hruid
    if fieldtype == 'integer':
        hru_ids = rootgrp.createVariable(hrufieldname,'i4',(hrufieldname))                # Coordinate Variable (32-bit signed integer)
    if fieldtype == 'real':
        hru_ids = rootgrp.createVariable(hrufieldname,'i4',(hrufieldname))                # Coordinate Variable (32-bit signed integer)
    elif fieldtype == 'string':
        hru_ids = rootgrp.createVariable(hrufieldname,str,(hrufieldname))                 # Coordinate Variable (string type character)

    # Create fixed-length variables
    weights = rootgrp.createVariable('weight', 'f8', (hrufieldname,'overlapPoly'),fill_value=-999.)  # (64-bit floating point)
    overlapPolyIds = rootgrp.createVariable('overlapPolyId', 'i4', (hrufieldname,'overlapPoly'),fill_value=-999)
    lats = rootgrp.createVariable('latitude', 'f8',(hrufieldname,'overlapPoly'),fill_value=-999.)    # (64-bit floating point)
    lons = rootgrp.createVariable('longitude', 'f8',(hrufieldname,'overlapPoly'),fill_value=-999.)   # (64-bit floating point)
    overlapps = rootgrp.createVariable('overlaps','i4',(hrufieldname))               # 32-bit signed integer

    # Set variable attribute 
    # Set units
    lats.units = 'degrees north'
    lons.units = 'degrees east'
    # Set longname
    hru_ids.longname = 'HRU ID'
    overlapPolyIds.longname = 'Overlapping polygon ID'
    weights.longname = 'Areal weight of overlapping polygon for hru'
    lons.longname = 'Longitude of overlapping polygon centroid'
    lats.longname = 'Latitude of overlapping polygon centroid'
    overlapps.longname = 'Number of overlapping polygon'

    # Fill in global attributes
    rootgrp.history = 'Created %s' %time.ctime()

    # Writing data
    # "hruid" dimension
    # sort hru_ids
    keylist = in_weight.keys()
    keylist.sort()
    if fieldtype == 'integer':
        #hru_ids[:] = np.array(in_weight.keys())                         # Ths method works for int-type netcdf variable
        hru_ids[:] = np.array(keylist)                         # Ths method works for int-type netcdf variable
    elif fieldtype == 'string':
        #hru_ids[:] = np.array(in_weight.keys(), dtype=np.object)     # This method works for a string-type netcdf variable
        hru_ids[:] = np.array(keylist, dtype.np.object)                         # Ths method works for int-type netcdf variable
    #hru_ids[:] = np.array([str(x) for x in weights.keys()], dtype=np.object)   # This method will convert to string-type netcdf variable from int

    # Write variables in netCDF 
    # "overlaps" Variables
    numOverlapPoly_tmp = np.array([len(x) for x in in_weight.values()])
    #overlapps[:] = numOverlapPoly

    # "overlapPolyIds", "weight", "latitude" and "longitude" Variables
    #initilize array
    s = (len(hru_ids),max(numOverlapPoly_tmp))
    overlapPolyIdsArray = np.ones(s)*-999
    weightsArray        = np.ones(s)*-999.0
    lonsArray           = np.ones(s)*-999.0
    latsArray           = np.ones(s)*-999.0
    numOverlapPolyArray = np.ones(len(hru_ids))*-999

    #Fill array
    i=0
    for hruids in sorted(in_weight.iterkeys()):                                 # For each basin
        weightslist       = [weight for weight in in_weight[hruids]]           # Generate list of gridcell IDs and spatial weights
        j=0
        for x in weightslist:                                                   # For each grid cell in each basin
            overlapPolyIdsArray[i,j] = x[0]
            weightsArray[i,j] = x[1]
            lonsArray[i,j] = in_LatLon[x[0]][0]
            latsArray[i,j] = in_LatLon[x[0]][1]
            j = j+1
        numOverlapPolyArray[i] = j
        i = i+1
    #Write netCDf

    overlapPolyIds[:,:] = overlapPolyIdsArray
    weights[:,:]        = weightsArray
    lons[:,:]           = lonsArray
    lats[:,:]           = latsArray
    overlapps[:]        = numOverlapPolyArray

    # Close file
    rootgrp.close()

def main():
    # Variables
    pass

if __name__ == '__main__':
    print "Cannot run as __main__ !"
    #main()
