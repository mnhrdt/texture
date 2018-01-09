#!/binbash
# auteur : mariedautume

# exemple input bash scaledim2colouredurve.sh 01
set -e
set -x
im=$1
thresh=$2
. config.sh

bin/ncc_compute_shift \
    data/mcdsm/outdir_01_02/outdir_from_01_02/cdsm.tif\
    data/Challenge1_Lidar_nan.tif \
    > data/ncc_shift/ncc_shift_$im.txt

read x y s z < data/ncc_shift/ncc_shift_$im.txt

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

