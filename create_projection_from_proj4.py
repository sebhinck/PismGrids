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
    #print "##################"
    #print inSpatialRef.ExportToPrettyWkt() # display information on orig projection              
    #print "##################"
     
    x=np.arange(((xLim[0] + dx/2.)*1000.), (xLim[1]*1000.), (dx*1000.))
    y=np.arange(((yLim[0] + dx/2.)*1000.), (yLim[1]*1000.), (dx*1000.))

    xx = np.expand_dims(x,axis=0).repeat(len(y),axis=0).reshape(len(x)*len(y))
    yy = np.expand_dims(y,axis=0).repeat(len(x),axis=0).transpose().reshape(len(x)*len(y))  

    out=np.array(coordTransform.TransformPoints(np.array((xx,yy)).transpose())).reshape((len(x),len(y),3))
    outfile=nc.Dataset(filename, "w", format='NETCDF3_64BIT')                                   
    fdx = outfile.createDimension("x", len(x))                                                            
    fdy = outfile.createDimension("y", len(y))                                                            

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

    elif "Lambert_Azimuthal_Equal_Area" in name:
        name = 'lambert_azimuthal_equal_area'
        mapping.grid_mapping_name = name ;
        mapping.longitude_of_projection_origin = inSpatialRef.GetProjParm('longitude_of_center') ;
        mapping.latitude_of_projection_origin = inSpatialRef.GetProjParm('latitude_of_center') ;

    outfile.proj4 = proj4

    fvx = outfile.createVariable("x","f4",("x",))
    fvy = outfile.createVariable("y","f4",("y",))
    fvx[:] = x                                   
    fvy[:] = y                                   
    fvx.units = "m"
    fvx.axis = "X"
    fvy.units = "m"
    fvy.axis = "Y"

    lon = outfile.createVariable("lon","f4",("y","x"))
    lat = outfile.createVariable("lat","f4",("y","x"))
    ok_topg = outfile.createVariable("ocean_kill_topg","f4",("y","x"))
    ok_thk  = outfile.createVariable("ocean_kill_thk" ,"f4",("y","x"))

    lon.units = "degreesE"
    lon.long_name = "Longitude"
    lon.standard_name = "longitude"

    lat.units = "degreesN"
    lat.long_name = "Latitude"
    lat.standard_name = "latitude"

    ok_topg.units = "m"
    ok_topg.standard_name = "bedrock_surface_elevation"
    ok_topg.coordinates  = "lon lat"
    ok_topg.grid_mapping = "mapping"

    ok_thk.units = "m"
    ok_thk.standard_name = "land_ice_thickness"
    ok_thk.coordinates  = "lon lat"
    ok_thk.grid_mapping = "mapping"

    lon[:] = out[:,:,0]
    lat[:] = out[:,:,1]

    dx_margin = 1

    thk_arr = np.ones(ok_thk.shape) * 100.
    thk_arr[(np.arange(0,dx_margin), np.arange(-dx_margin, 0)), :] = 0.0
    thk_arr[:, (np.arange(0,dx_margin), np.arange(-dx_margin, 0))] = 0.0

    ok_topg[:] = -1
    ok_thk[:]  = thk_arr[:,:]

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
