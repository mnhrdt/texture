#!/binbash
# auteur : mariedautume

# exemple input bash scaledim2colouredurve.sh 01
set -e
set -x
im=$1

. config.sh

bin/ncc_compute_shift \
    data/mcdsm/outdir_$im\_*/outdir_from_$im\_*/cdsm.tif\
    data/Challenge1_Lidar_nan.tif \
    | read x y s z
#     < data/ncc/ncc_shift_$im.txt
#
# read x y s z < data/ncc/ncc_shift_$im.txt

bin/zbuffer \
    data/mesh_curve_scaled_remeshed_02.off\
    $ZONE $ORIG_X $ORIG_Y $SCALE_X $SCALE_Y\
    data/PAN/img_$im.rpc \
    data/xywh/xywh_$im.txt \
    exp/soutput/matches/matches_$im.tif \
    -ox $x -oy $y -oz $z -xmin $CROP_X -ymin $CROP_Y --res $RES_MESH

bin/colormultiple_mesh \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    data/PAN/img_$im.ntf \
    data/PAN/img_$im.rpc \
    data/MSI/img_$im.ntf \
    data/MSI/img_$im.xml \
    exp/soutput/matches/matches_$im.tif  \
    exp/soutput/vc_$im.tif

qauto -i -f -p 0.5 exp/soutput/vc_$im.tif exp/soutput/vc_$im.tiff
# qauto -i -f exp/soutput/vc_$im.tif exp/soutput/vc_$im.tiff

bin/write_coloured_ply \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    exp/soutput/vc_$im.tiff \
    exp/soutput/new_ccoloured_$im\_02.ply 

