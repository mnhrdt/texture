#!/binbash
# auteur : mariedautume

# exemple input bash all_new.sh 01
set -e
set -x
im=$1

. scripts/config.sh

read x y s z < data/ncc_shift/ncc_shift_$im.txt
#x=`sed -n '1p' < data/martin/martin_$im.txt | cut -f 2 -d \( | cut -f 1 -d ,`
#y=`sed -n '1p' < data/martin/martin_$im.txt | cut -f 2 -d , | cut -f 1 -d \)`
#z=`sed -n '3p' < data/martin/martin_$im.txt`

#echo $x 
#echo $y
#echo $z


mkdir -p output
az=`gdalinfo data/MSI/img_$im.ntf | grep SUN_AZI | cut -f 2 -d=`
el=`gdalinfo data/MSI/img_$im.ntf | grep SUN_ELE | cut -f 2 -d=`

#bin/triproc data/mesh_curve_scaled_remeshed_02.off output/edges.txt

bin/get_utm_normal_shadow \
    data/mesh_curve_scaled_remeshed_02.off \
    $SCALE_X $SCALE_Y \
    $ORIG_X $ORIG_Y \
    $ZONE \
    data/PAN/img_$im.rpc \
    $az $el \
    output/utm_coord_$im.tif \
    output/theoric_sun_$im.tif \
    output/scalars_$im.tif \
    -ox $x -oy $y -oz $z -xmin $CROP_X -ymin $CROP_Y 

bin/colorize_vertices_from_one_image \
    data/mesh_curve_scaled_remeshed_02.off \
    data/PAN/img_$im.ntf \
    data/PAN/img_$im.rpc \
    data/MSI/img_$im.ntf \
    data/MSI/img_$im.xml \
    -21 \
    output/utm_coord_$im.tif \
    output/theoric_sun_$im.tif \
    output/pan_$im.tif \
    output/msi_$im.tif \
    output/rgb_$im.tif \
    output/real_sun_$im.tif 


#qauto -i -f -p 0.5 exp/soutput/vc_$im.tif exp/soutput/vc_$im.tiff
# qauto -i -f exp/soutput/vc_$im.tif exp/soutput/vc_$im.tiff

#bin/write_coloured_ply \
#    exp/output/mesh_curve_scaled_remeshed_02.off \
#    exp/soutput/vc_$im.tiff \
#    exp/soutput/mesh_registred_ps_$im\_02.ply 

