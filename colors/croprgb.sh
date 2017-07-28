#!/bin/bash

# This is a script for extracting a RGB crop from a whole PAN+MSI image.
# INPUT: 4 numbers of the crop
#        the names PAN and MSI files
#        the offset correction of the MSI (in MSI pixel units (!) )
# OUTPUT: a 8-bit rgb image
#
# required binaries: gdal_translate, plambda, upsa
#
# caveat: right now it only works for Worldview3 images (8-band)

# set debugging flags
set -e
set -x

# check number of input arguments
if test $# -ne 9 ; then
	printf "usage:\n\tcroprgb.sh x0 y0 w h ox oy in_pan in_msi out_rgb\n"
	#                            1  2  3 4 5  6  7      8      9
	exit 1
fi

# assign input arguments to variables with reasonable names
X=$1
Y=$2
W=$3
H=$4
OX=$5
OY=$6
PAN=$7
MSI=$8
RGB=$9

# setup the crop strings for each image
CROP_PAN="$X $Y $W $H"
C1=`echo "$OX+$X/4" | bc`
C2=`echo "$OY+$Y/4" | bc`
C3=`echo "1+$W/4" | bc`
C4=`echo "1+$H/4" | bc`
CROP_MSI="$C1 $C2 $C3 $C4"

# create temporary directory for intermediate files
TPD=`mktemp -d`
echo "crop pan $CROP_PAN crop msi $CROP_MSI"
# extract each crop using gdal_translate
gdal_translate -ot uint16 -srcwin $CROP_PAN $PAN $TPD/pan.tif
gdal_translate -ot uint16 -srcwin $CROP_MSI $MSI $TPD/msi.tif
gdal_translate -ot uint16 -srcwin $CROP_PAN $PAN pan.tif
gdal_translate -ot uint16 -srcwin $CROP_MSI $MSI msi.tif

# obtain rgb from msi (approximate formula)
MSITORGB="x[4] x[2] 0.8 * x[5] 0.1 * + x[1] 1.2 * join3 log 5 7.3 range"
plambda $TPD/msi.tif "$MSITORGB" | upsa 4 1 - $TPD/rgb.tif

# compute P+XS (already in rgb space)
plambda $TPD/pan.tif $TPD/rgb.tif "dup vnorm / * -30 600 qe" -o $RGB

# delete temporary directory
rm -rf $TPD
