#!/bin/bash
# auteur : mariedautume

# normalement appelé depuis entire.sh
# sinon bash dir config IM CDSM PAN_ntf PAN_rpc MSI_ntf PSI_rpc
set -e
#set -x

dir=$1         # folder where to print the results
config=$2      # configuration file (exemple scripts/config.sh
IM=$3          # image id typically a number
CDSM=$4        # CDSM corresponding to the image 
PAN_ntf=$5     # satellite image (panchromatic)  
PAN_rpc=$6     # panchromatic rpc                
MSI_ntf=$7     # satellite image (MSI)           
MSI_rpc=$8     # MSI rpc                         

# read configaration file
. $config

# creating temporary files if they don't exist
mkdir -p $dir/ncc_shift
mkdir -p $dir/utm_coord 
mkdir -p $dir/theoric_sun
mkdir -p $dir/real_sun
mkdir -p $dir/scalars
mkdir -p $dir/pan
mkdir -p $dir/msi
mkdir -p $dir/rgb

# computing shift with reference dsm
#$RECALAGE $CDSM $DSM $dir/ncc_shift/shifted_$IM.tif > $dir/ncc_shift/ncc_shift_$IM.txt
$RECALAGE $CDSM $DSM $dir/ncc_shift/shifted_$IM.tif > $dir/ncc_shift/ncc_shift_$IM.txt
read x y z < $dir/ncc_shift/ncc_shift_$IM.txt

# get azimuth et elevation from the msi images
az=`gdalinfo $MSI_ntf | grep SUN_AZI | cut -f 2 -d=`
el=`gdalinfo $MSI_ntf | grep SUN_ELE | cut -f 2 -d=`

# get for each vertex its utm coordinates if visible (nan otherwise)
#                         theoric shadow from the geometry
#                         scalar product between satellite orientation and
#                                surface normal
bin/get_utm_normal_shadow \
    $MESH \
    $SCALE_X $SCALE_Y \
    $ORIG_X $ORIG_Y \
    $ZONE \
    $PAN_rpc \
    $az $el \
    $dir/utm_coord/utm_coord_$IM.tif \
    $dir/theoric_sun/theoric_sun_$IM.tif \
    $dir/scalars/scalars_$IM.tif \
    -ox $x -oy $y -oz $z -xmin $CROP_X -ymin $CROP_Y 

# get for each vertex its pan, msi and rgb values and shadow obtained from 
# thresholding
bin/colorize_vertices_from_one_image \
    $MESH \
    $PAN_ntf \
    $PAN_rpc \
    $MSI_ntf \
    $MSI_rpc \
    $ZONE \
    $dir/utm_coord/utm_coord_$IM.tif \
    $dir/theoric_sun/theoric_sun_$IM.tif \
    $dir/pan/pan_$IM.tif \
    $dir/msi/msi_$IM.tif \
    $dir/rgb/rgb_$IM.tif \
    $dir/real_sun/real_sun_$IM.tif 

# quantification to be able to directly look at the intermediary results 
# with mflip
########### TO DO ##############
# better choice of quantization
################################

qauto -i -f -p 0.5 $dir/rgb/rgb_$IM.tif $dir/rgb/rgb_$IM.tiff

