#!/binbash
# auteur : mariedautume

# exemple input bash im2colouredurve.sh 01
set -e
set -x
im=$1
thresh=$2

read x y s z < data/ncc_shift/ncc_shift_$im.txt

az=`gdalinfo data/MSI/img_$im.ntf | grep SUN_AZI | cut -f 2 -d=`
el=`gdalinfo data/MSI/img_$im.ntf | grep SUN_ELE | cut -f 2 -d=`

bin/get_corners_utm data/lidar_curve.tif -21 $az $el > corners.txt

bin/zbuffer \
    data/scale_$im.txt \
    -21 \
    data/PAN/img_$im.rpc \
    corners.txt \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    exp/soutput/vs_$im.tif \
    $az $el -s 1 \
    --xywh data/xywh/xywh_$im.txt \
    --proj exp/output/shadow_proj_$im.png \
    -ox $x -oy $y -oz $z -xmin 1505 -ymin 891 --res 0.45

bin/write_coloured_ply \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    exp/soutput/vs_$im.tif \
    exp/soutput/shadow_$im\_02.ply 

