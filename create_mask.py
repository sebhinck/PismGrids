#!/usr/bin/env python                                                                                       

import osr
import netCDF4 as nc
import numpy as np  
import sys          
import getopt       
from argparse import ArgumentParser


def create_projection (filename, proj4, xLim, yLim, dx, opts={}):
    inSpatialRef = osr.SpatialReference() # orig projection
    inSpatialRef.ImportFromProj4(proj4)
    outSpatialRef = osr.SpatialReference() # target projection
    outSpatialRef.ImportFromEPSG(4326) # WGS 1984
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef) # create converter
    print "##################"
    print inSpatialRef.ExportToPrettyWkt() # display information on orig projection              
    print "##################"
     
    x=np.arange(((xLim[0] + dx/2.)*1000.), (xLim[1]*1000.), (dx*1000.))
    y=np.arange(((yLim[0] + dx/2.)*1000.), (yLim[1]*1000.), (dx*1000.))

    xx = np.expand_dims(x,axis=0).repeat(len(y),axis=0).reshape(len(x)*len(y)) #.transpose()
    yy = np.expand_dims(y,axis=0).repeat(len(x),axis=0).transpose().reshape(len(x)*len(y))  

    out=np.array(coordTransform.TransformPoints(np.array((xx,yy)).transpose())).reshape((len(x),len(y),3))
    outfile=nc.Dataset(filename, "w", format='NETCDF3_64BIT')                                   
    fdx = outfile.createDimension("x", len(x))                                                            
    fdy = outfile.createDimension("y", len(y))                                                            
    fdy = outfile.createDimension("grid_corners", 4)      

    mapping = outfile.createVariable("mapping","f4")               
    mapping.ellipsoid = inSpatialRef.GetAttrValue('GEOGCS').replace(" ", "") ;
    mapping.false_easting = inSpatialRef.GetProjParm('false_easting') ;
    mapping.false_northing = inSpatialRef.GetProjParm('false_northing') ;
    
    name = inSpatialRef.GetAttrValue('PROJECTION')
    if "Lambert_Conformal_Conic" in name:
        name = 'lambert_conformal_conic'
        mapping.grid_mapping_name = name ;
        mapping.latitude_of_projection_origin = inSpatialRef.GetProjParm('latitude_of_origin') ;
        mapping.longitude_of_central_meridian = inSpatialRef.GetProjParm('central_meridian') ;
        mapping.standard_parallel = [inSpatialRef.GetProjParm('standard_parallel_1'), inSpatialRef.GetProjParm('standard_parallel_2')]
        
    elif "Stereographic" in name:
        name = 'polar_stereographic'
        mapping.grid_mapping_name = name ;
	lat0 = inSpatialRef.GetProjParm('latitude_of_origin') ;
	#Dirty fix...
	if lat0 > 0:
		lat0 = 90;	
	else:
		lat0 = -90;	
        mapping.latitude_of_projection_origin = lat0 ;
        mapping.standard_parallel = inSpatialRef.GetProjParm('latitude_of_origin') ;
        mapping.straight_vertical_longitude_from_pole = inSpatialRef.GetProjParm('central_meridian') ;

    outfile.proj4 = proj4

    fvx = outfile.createVariable("x","f4",("x",))
    fvy = outfile.createVariable("y","f4",("y",))
    fvx[:] = x                                   
    fvy[:] = y                                   
    fvx.units = "m"                              
    fvy.units = "m"   

    lon = outfile.createVariable("lon","f4",("y","x"), fill_value=-9.e9)
    lat = outfile.createVariable("lat","f4",("y","x"), fill_value=-9.e9)

    lon.units = "degrees"
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon.bounds = "grid_corner_lon" 
    lon._CoordinateAxisType = "Lon"
    lon.grid_mapping = "mapping"

    lat.units = "degrees"
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat.bounds = "grid_corner_lat"
    lat._CoordinateAxisType = "Lat"
    lat.grid_mapping = "mapping"

    lon[:] = out[:,:,0]
    lat[:] = out[:,:,1]
    xoff=(x[1]-x[0])/2.
    yoff=(y[1]-y[0])/2.
    addx = [ xoff, xoff,-xoff,-xoff]
    addy = [-yoff, yoff, yoff,-yoff]

    grid_corner_lat = outfile.createVariable("grid_corner_lat","f4",("y", "x", "grid_corners"), fill_value=-9.e9)
    grid_corner_lon = outfile.createVariable("grid_corner_lon","f4",("y", "x", "grid_corners"), fill_value=-9.e9)

    grid_corner_lat.units = "degrees"
    grid_corner_lon.units = "degrees"

    for i in xrange(4):
        corner_temp = np.array(coordTransform.TransformPoints(np.array((xx + addx[i], yy + addy[i])).transpose())).reshape((len(x),len(y),3))
        grid_corner_lon[:,:,i] = corner_temp[:, :, 0]
        grid_corner_lat[:,:,i] = corner_temp[:, :, 1]

    outfile.close()
#################


def parse_args():
  parser = ArgumentParser()
  parser.description = "compare slopes of two variables from two files"
  parser.add_argument("FILE")
  parser.add_argument("proj4")
  parser.add_argument("x_lim", nargs=2, type=float)
  parser.add_argument("y_lim", nargs=2, type=float)
  parser.add_argument("dx", type=float)
  parser.add_argument("-v", "--verbose",
                    help='''Be verbose''', action="store_true")

  options = parser.parse_args()
  return options


def main():
  options = parse_args()
  if options.verbose:
    print (dir(options))
    print options.FILE

  create_projection(options.FILE, options.proj4 ,options.x_lim, options.y_lim, options.dx, vars(options))



if __name__ == "__main__":
    main()
