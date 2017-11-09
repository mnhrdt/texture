#!/binbash
# auteur : mariedautume

# exemple bash shadow.sh 01
set -e
set -x
im=$1

az=`gdalinfo data/MSI/img_$im.ntf | grep SUN_AZI | cut -f 2 -d=`
el=`gdalinfo data/MSI/img_$im.ntf | grep SUN_ELE | cut -f 2 -d=`

bin/get_corners_utm \
    data/lidar_curve.tif \
    -21 \
    $az $el \
    > corners.txt

bin/shadow \
     data/lidar_curve.tif \
     -21 \
     $az $el \
     corners.txt \
     data/mesh_curve_remeshed.off \
     vs.tif \
     -res 0.45 

bin/write_coloured_ply \
    data/mesh_curve_remeshed.off \
    vs.tif  \
    shadow_$im.ply


