#!/usr/local/bin/bash
# auteur : mariedautume

# exemple input bash all_new.sh 01
set -e
set -x

dir=$1
config=$2
IM=$3
CDSM=$4
PAN_ntf=$5
PAN_rpc=$6
MSI_ntf=$7
MSI_rpc=$8

. $config

$RECALAGE $CDSM $DSM > $dir/ncc_shift/ncc_shift_$IM.txt
read x y z < $dir/ncc_shift/ncc_shift.txt

mkdir -p $dir/utm_coord 
mkdir -p $dir/theoric_sun
mkdir -p $dir/real_sun
mkdir -p $dir/scalars
mkdir -p $dir/pan
mkdir -p $dir/msi
mkdir -p $dir/rgb

az=`gdalinfo $MSI_ntf | grep SUN_AZI | cut -f 2 -d=`
el=`gdalinfo $MSI_ntf | grep SUN_ELE | cut -f 2 -d=`

#bin/triproc data/mesh_curve_scaled_remeshed_02.off output/edges.txt
    #tmp/mesh/remeshed_dsm_37.off \

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

    #tmp/mesh/remeshed_dsm_37.off \
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


qauto -i -f -p 0.5 $dir/rgb/rgb_$IM.tif $dir/rgb/rgb_$IM.tiff
# qauto -i -f exp/soutput/vc_$im.tif exp/soutput/vc_$im.tiff

#bin/write_coloured_ply \
#    exp/output/mesh_curve_scaled_remeshed_02.off \
#    exp/soutput/vc_$im.tiff \
#    exp/soutput/mesh_registred_ps_$im\_02.ply 

