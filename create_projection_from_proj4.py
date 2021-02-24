#!/usr/bin/env python

import osr
import netCDF4 as nc
import xarray as xr
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

    x=np.arange(((xLim[0] + dx/2.)*1000.), (xLim[1]*1000.), (dx*1000.))
    y=np.arange(((yLim[0] + dx/2.)*1000.), (yLim[1]*1000.), (dx*1000.))

    xa = xr.DataArray(data=np.zeros([len(x), len(y)]), coords = [('x',x), ('y', y)])

    def ctransform(a):
        yi, xi = xr.broadcast(a.y, a.x)

        def f(x,y,z=0):
            res = coordTransform.TransformPoint(x,y,z)

            return np.stack(res, axis=-1)

        res = xr.apply_ufunc(f, xi, yi, vectorize=True, output_core_dims= [["q"]])

        lon = res.isel(q=0)
        lat = res.isel(q=1)

        return lon, lat

    lon, lat = ctransform(xa)

    lon.attrs['units'] = "degreesE"
    lon.attrs['long_name'] = "Longitude"
    lon.attrs['standard_name'] = "longitude"

    lat.attrs['units'] = "degreesN"
    lat.attrs['long_name'] = "Latitude"
    lat.attrs['standard_name'] = "latitude"

    ds = xr.Dataset()

    ds['lon'] = lon
    ds['lat'] = lat

    ds['mapping'] = xr.DataArray(None)
    ds['mapping'].attrs = dict(ellipsoid = inSpatialRef.GetAttrValue('GEOGCS').replace(" ", ""),
                               false_easting = inSpatialRef.GetProjParm('false_easting'),
                               false_northing = inSpatialRef.GetProjParm('false_northing'))

    name = inSpatialRef.GetAttrValue('PROJECTION')
    if "Lambert_Conformal_Conic" in name:
        name = 'lambert_conformal_conic'
        ds['mapping'].attrs['grid_mapping_name'] = name ;
        ds['mapping'].attrs['latitude_of_projection_origin'] = inSpatialRef.GetProjParm('latitude_of_origin') ;
        ds['mapping'].attrs['longitude_of_central_meridian'] = inSpatialRef.GetProjParm('central_meridian') ;
        ds['mapping'].attrs['standard_parallel'] = [inSpatialRef.GetProjParm('standard_parallel_1'), inSpatialRef.GetProjParm('standard_parallel_2')]

    elif "Stereographic" in name:
        name = 'polar_stereographic'
        ds['mapping'].attrs['grid_mapping_name'] = name ;
        lat0 = inSpatialRef.GetProjParm('latitude_of_origin') ;
        #Dirty fix...
        if lat0 > 0:
          lat0 = 90;
        else:
          lat0 = -90;
        ds['mapping'].attrs['latitude_of_projection_origin'] = lat0 ;
        ds['mapping'].attrs['standard_parallel'] = inSpatialRef.GetProjParm('latitude_of_origin') ;
        ds['mapping'].attrs['straight_vertical_longitude_from_pole'] = inSpatialRef.GetProjParm('central_meridian') ;

    elif "Lambert_Azimuthal_Equal_Area" in name:
        name = 'lambert_azimuthal_equal_area'
        ds['mapping'].attrs['grid_mapping_name'] = name ;
        ds['mapping'].attrs['longitude_of_projection_origin'] = inSpatialRef.GetProjParm('longitude_of_center') ;
        ds['mapping'].attrs['latitude_of_projection_origin'] = inSpatialRef.GetProjParm('latitude_of_center') ;

    dx_margin = 1

    ds['ocean_kill_topg'] = xr.DataArray(np.ones_like(lon) * -1,
                                         dims=lon.dims,
                                         attrs=dict(
                                           units = 'm',
                                           standard_name = "bedrock_altitude",
                                           coordinates  = "lon lat",
                                           grid_mapping = "mapping"
                                         )
                                        )

    tmp = np.ones_like(lon) * 100.
    tmp[(np.arange(0,dx_margin), np.arange(-dx_margin, 0)), :] = 0.0
    tmp[:, (np.arange(0,dx_margin), np.arange(-dx_margin, 0))] = 0.0
    ds['ocean_kill_thk'] = xr.DataArray(tmp, 
                                        dims=lon.dims,
                                        attrs=dict(
                                          units = 'm',
                                          standard_name = "land_ice_thickness",
                                          coordinates  = "lon lat",
                                          grid_mapping = "mapping"
                                        )
                                       )

    tmp = np.ones_like(lon, dtype=int)
    tmp[(np.arange(0,dx_margin), np.arange(-dx_margin, 0)), :] = 0
    tmp[:, (np.arange(0,dx_margin), np.arange(-dx_margin, 0))] = 0
    ds['land_ice_area_fraction_retreat'] = xr.DataArray(tmp,
                                                        dims=lon.dims,
                                                        attrs=dict(
                                                        units = '1',
                                                        standard_name = "land_ice_area_fraction_retreat",
                                                        coordinates  = "lon lat",
                                                        grid_mapping = "mapping"
                                                       )
                                                      )
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
