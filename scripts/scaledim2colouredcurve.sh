#!/binbash
# auteur : mariedautume

# exemple input bash im2colouredurve.sh 01
set -e
set -x
im=$1
thresh=$2

read x y s z < data/ncc_shift/ncc_shift_$im.txt

bin/zbuffer \
    data/scale_$im.txt \
    -21 \
    data/PAN/img_$im.rpc \
    data/xywh/xywh_$im.txt \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    exp/soutput/matches/matches_$im.tif \
    -ox $x -oy $y -oz $z -xmin 1505 -ymin 891 --res 0.45

bin/colormultiple_mesh \
    exp/output/mesh_curve_scaled_remeshed_02.off \
    data/xywh/xywh_$im.txt \
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

