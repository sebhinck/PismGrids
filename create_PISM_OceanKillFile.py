#!/usr/bin/env python                                                                                       

import osr
import netCDF4 as nc
import numpy as np  
import sys          
import getopt       
from argparse import ArgumentParser


def   create_PISM_OceanKillFile(in_file, dx ,out_file):

    src = nc.Dataset(in_file,  "r")
    dst = nc.Dataset(out_file, "w")

    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
       x = dst.createVariable(name, variable.datatype, variable.dimensions)
       dst[name][:] = src[name][:]
       # copy variable attributes all at once via dictionary
       dst[name].setncatts(src[name].__dict__)

    src.close()

    topg = dst.createVariable("topg","f4",("y","x"), fill_value=-9.e9)
    topg.units = "m"
    topg.coordinates = "lat lon"
    topg.grid_mapping = "mapping"
    topg.standard_name = "bedrock_surface_elevation"

    thk = dst.createVariable("thk","f4",("y","x"), fill_value=-9.e9)
    thk.units = "m"
    thk.coordinates = "lat lon"
    thk.grid_mapping = "mapping"
    thk.standard_name = "land_ice_thickness"

    thk_arr = np.ones(thk.shape) * 100.
    thk_arr[(np.arange(0,dx), np.arange(-dx, 0)),:] = 0.0
    thk_arr[:, (np.arange(0,dx), np.arange(-dx, 0))] = 0.0

    thk[:]  = thk_arr[:,:]
    topg[:] = -1.

    dst.close()

#################


def parse_args():
  parser = ArgumentParser()
  parser.description = "create PISM ocean kill file"
  parser.add_argument("in_file")
  parser.add_argument("dx", type=int)
  parser.add_argument("out_file")
  parser.add_argument("-v", "--verbose",
                    help='''Be verbose''', action="store_true")

  options = parser.parse_args()
  return options


def main():
  options = parse_args()
  if options.verbose:
    print (dir(options))
    print options.in_file
    print options.dx
    print options.out_file

  create_PISM_OceanKillFile(options.in_file, options.dx ,options.out_file)



if __name__ == "__main__":
    main()
