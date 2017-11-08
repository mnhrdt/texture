#!/binbash
# auteur : mariedautume

# exemple input bash im2curve.sh 01
set -e
set -x
im=$1
thresh=$2

read x y s z < data/ncc_shift/ncc_shift_$im.txt

bin/zbuffer \
    data/mcdsm/outdir_$im*/mcdsm.tif \
    -21 \
    exp/data/crop/crop_img_$im.tif  \
    data/PAN/img_$im.rpc \
    exp/output/xywh/xywh_$im.txt \
    data/mesh_curve_remeshed.off \
    exp/output/matches/matches_$im.tif \
    -ox $x -oy $y -oz $z -xmin 1505 -ymin 891 --res 0.45

bin/colormultiple_mesh \
    exp/output/mesh_curve_remeshed.off \
    exp/data/crop/crop_img_$im.tif \
    data/MSI/img_$im.ntf \
    data/MSI/img_$im.xml \
    exp/output/matches/matches_$im.tif  \
    exp/output/vc_$im.tif

qauto exp/output/vc_$im.tif exp/output/vc_$im.tiff

bin/write_coloured_ply \
    exp/output/mesh_curve_remeshed.off \
    exp/output/vc_$im.tif \
    exp/output/ccoloured_$im.ply \
    -t $thresh

