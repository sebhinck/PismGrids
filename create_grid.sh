#!/bin/bash

#run:
# ./generate_mask.sh dx_km domain_name left_km right_km bottom_km top_km proj4_string out_folder

#This describes the projection used for the PISM Laurentide domain
default_proj4="+proj=lcc +lat_1=45 +lat_2=65 +lat_0=63 +lon_0=-100 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

dx=${1:-"20"}
dom=${2:-"laurentide"}
x0=${3:-"-3900"}
x1=${4:-"3900"}
y0=${5:-"-3300"}
y1=${6:-"3300"}
proj4=${7:-$default_proj4}
outfolder=${8:-"/home/ollie/shinck/projects/PISM/grids"}

outfolder=${outfolder}/${dom}/
mkdir -p $outfolder

name=pismr_${dom}_${dx}km
gfile=${name}.nc

rm -f $gfile
./create_projection_from_proj4.py ${gfile} "${proj4}" $x0 $x1 $y0 $y1 $dx

#PISM Bin dir need to be set in PATH
nc2cdo.py $gfile

#generate griddes
cdo griddes $gfile > ${name}.griddes

#extract x,y
ncks -O -v x,y,lon,lat,mapping $gfile ${name}_xy.nc

echo $(ncks -M -m $gfile | grep -E -i "^global attribute [0-9]+: proj4" | cut -f 11- -d ' ') > ${name}.proj4

border_width=2
./create_PISM_OceanKillFile.py ${name}_xy.nc  $border_width  ${name}_ocean_kill_file.nc

cp ${name}.griddes $outfolder
cp ${name}_xy.nc $outfolder
cp ${name}.proj4 $outfolder
cp ${name}_ocean_kill_file.nc $outfolder

rm ${name}.griddes ${name}_xy.nc ${name}.proj4 ${name}.nc ${name}_ocean_kill_file.nc

