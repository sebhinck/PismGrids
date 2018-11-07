#!/bin/bash

dom=$1
dx=$2
outfolder=${3:-"/home/ollie/shinck/projects/PISM/grids"}

module purge
module load pism/dev_intel_impi netcdf-tools cdo gdal python

case "$dom" in
  nhem)
    #EPSG3413
    proj4="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    xmin=-6000
    xmax=6000
    ymin=-6000
    ymax=6000
  ;;
  greenland)
    #EPSG3413 - InitMIP 20 km, 10 km, 5 km, 2 km or 1 km
    #See: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_1_.E2.80.93_Output_grid_definition_and_interpolation
    proj4="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    xmin_center=-720
    xmax_center=960
    ymin_center=-3450
    ymax_center=-570

    xmin=$(echo "$xmin_center - (0.5 * $dx)" | bc)
    xmax=$(echo "$xmax_center + (0.5 * $dx)" | bc)
    ymin=$(echo "$ymin_center - (0.5 * $dx)" | bc)
    ymax=$(echo "$ymax_center + (0.5 * $dx)" | bc)
  ;;
  greenland_handbook)
    #Bamber
    #See: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Greenland#Appendix_1_.E2.80.93_Output_grid_definition_and_interpolation
    proj4="+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    xmin_center=-800
    xmax_center=700
    ymin_center=-3400
    ymax_center=-600

    xmin=$(echo "$xmin_center - (0.5 * $dx)" | bc)
    xmax=$(echo "$xmax_center + (0.5 * $dx)" | bc)
    ymin=$(echo "$ymin_center - (0.5 * $dx)" | bc)
    ymax=$(echo "$ymax_center + (0.5 * $dx)" | bc)
  ;;
  laurentide)
    proj4="+proj=lcc +lat_1=45 +lat_2=65 +lat_0=63 +lon_0=-100 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
    xmin=-3900
    xmax=3900
    ymin=-3300
    ymax=3300
  ;;
  antarctica)
    #ISMIP6 32km, 16 km, 8 km, 4 km, 2 km or 1 km
    #See: http://www.climate-cryosphere.org/wiki/index.php?title=InitMIP-Antarctica#Appendix_1_.E2.80.93_Output_grid_definition_and_interpolation
    proj4="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    xmin_center=-3040
    xmax_center=3040
    ymin_center=-3040
    ymax_center=3040

    xmin=$(echo "$xmin_center - (0.5 * $dx)" | bc)
    xmax=$(echo "$xmax_center + (0.5 * $dx)" | bc)
    ymin=$(echo "$ymin_center - (0.5 * $dx)" | bc)
    ymax=$(echo "$ymax_center + (0.5 * $dx)" | bc)
  ;;
  LIS_Evan)
    #The false easting and northing and the ranges of x and y correspond to the corner points at 5km resolution:
    #   west_latitude=25, west_longitude=-135, east_latitude=85, east_longitude=3
    proj4="+proj=laea +lon_0=-94 +lat_0=60 +units=m +x_0=4104009.407173712 +y_0=2625682.633840935"
    
    #At 5km resolution the cell centers of the corner cells are at
    #xmin_center=0
    #xmax_center=7750
    #ymin_center=0
    #ymax_center=5950

    xmin=-2.5
    xmax=7752.5
    ymin=-2.5
    ymax=5952.5
  ;;
  *)
    echo "Domain \"${dom}\" not defined!"
    exit 1
esac

echo
echo "#####################################"
echo "Preparing grid for domain \"${dom}\" at ${dx}km resolution"

Lx=$(echo "($xmax - $xmin)" | bc)
Ly=$(echo "($ymax - $ymin)" | bc)
Nx=$(echo "($Lx/$dx)" | bc)
Ny=$(echo "($Ly/$dx)" | bc)

echo "-> Grid: (${Lx}km x ${Ly}km) = (${Nx} x ${Ny}) @ ${dx}km resolution"

func_call="./create_grid.sh $dx $dom $xmin $xmax $ymin $ymax \"$proj4\" $outfolder"
echo $func_call
echo "#####################################"
echo

eval $func_call

