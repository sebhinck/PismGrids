#!/bin/sh

ifile=$1
ofile=$2

ncap2 -s mask=lon*0 $ifile temp
ncatted -a _CoordinateAxisType,mask,d,,, temp
ncatted -a coordinates,mask,o,c,'lat lon' temp
ncatted -a long_name,mask,d,,, temp
ncatted -a bounds,mask,d,,, temp
ncatted -a standard_name,mask,d,,, temp
ncatted -a bounds,lon,d,,, temp
ncatted -a bounds,lat,d,,, temp
ncks -v mask,lat,lon,x,y,mapping temp $ofile
rm temp

