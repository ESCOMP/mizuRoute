#!/usr/bin/env python
"""
 Name:        map_poly2poly.py
 Purpose:
 Author:       Kevin Sampson
               Associate Scientist
               National Center for Atmospheric Research
               ksampson@ucar.edu
 Created:      10/23/2013
 last Modified:02/02/2014 naoki

 Modified:01/14/2015 naoki
  added hru_id name for overlapping polygon for option
"""

# Import modules
import time
import getopt  
import os
import sys
import glob
import GDAL_OGR_Functions
tic = time.time()

def gridtobasin_function(infeatures, hrufieldname1, inpoly, hrufieldname2):

    # Input directories
    indir = os.path.dirname(__file__)                                               # Store current working directory
    
    # Intermediate shapefile
    outputFile1 = os.path.join(indir, os.path.basename(inpoly)[:-4]+"_proj.shp")
    # Final output netCDF 
    outputFile3 = os.path.join(indir, os.path.basename(infeatures)[:-4]+'.nc')

    # Delete projected shapefile if exist
    print "Removing existing file: %s" %outputFile1
    [os.remove(infile) for infile in glob.glob(os.path.join(indir, os.path.basename(inpoly)[:-4]+'_proj.*')) if os.path.isfile(infile)==True]
    if os.path.isfile(outputFile3) == True:
        print "Removing existing file: %s" %outputFile3
        os.remove(outputFile3)

    # Step 1: Compute centroid lat/lon for inpoly
    tic1 = time.time()
    latlondict=GDAL_OGR_Functions.polyXY(inpoly,hrufieldname2)
    print "Step 1 completed in %s seconds." %(time.time()-tic1)

    # Step 2: Project gridded polygons to coordinate system of input features
    tic1 = time.time()
    outputFile1 = GDAL_OGR_Functions.project_to_input(infeatures, outputFile1, inpoly)
    print "Step 2 completed in %s seconds." %(time.time()-tic1)

    # Step 3: Calculate per gridcell feature coverage fractions
    tic1 = time.time()
    spatialweights, fieldtype = GDAL_OGR_Functions.spatial_weights(infeatures,  hrufieldname1, outputFile1, hrufieldname2)
    print "Step 3 completed in %s seconds." %(time.time()-tic1)

    # Step 4: Write weights to ragged array netcdf file
    tic1 = time.time()
    GDAL_OGR_Functions.create_ncfile(latlondict, spatialweights, outputFile3, hrufieldname1, fieldtype)
    print "Step 4 completed in %s seconds." %(time.time()-tic1)

    print "Total time elapsed: %s seconds." %(time.time()-tic)

use = '''
Usage: %s <hruPoly.sh> <FieldName> <overlayPoly.sh>
      <hruPoly.sh>     -> hru Polygon Shapefile  
      <FieldName1>     -> Polygon identifier for hruPoly
      <overlayPoly.sh> -> overlaying polygon shapefile
      <FieldName2>     -> Polygon identifier for overlayPoly
'''
if __name__ == '__main__':

  def usage():
      sys.stderr.write(use % sys.argv[0])
      sys.exit(1)
  try:
      (opts, args) = getopt.getopt(sys.argv[1:], 'h')
  except getopt.error:
      usage()

  if len(sys.argv) != 5:
      usage()
  else: 

    # Input arguments
    hruPoly = sys.argv[1]
    hrufieldname1 = sys.argv[2] # usually "hru_id" 
    overlayPoly   = sys.argv[3]
    hrufieldname2 = sys.argv[4] # usually "hru_id" 
      
    gridtobasin_function(hruPoly, hrufieldname1, overlayPoly, hrufieldname2)
   
    # Clean up the files
    [os.remove(rmfile) for rmfile in glob.glob(os.path.basename(overlayPoly)[:-4]+'_proj.*')]

