#!/binbash
# auteur : mariedautume

# exemple input bash meshmoche.sh /home/coco/src/s2p-iarpa/input/ essai_small/ data/small_lidar.tif
set -e
set -x
im=$1

bin/get_corners data/mcdsm/outdir_$im*/mcdsm.tif -21 data/PAN/img_$im.rpc  > tmp/output/xywh/xywh_$im.txt

read x y w h < tmp/output/xywh/xywh_$im.txt

gdal_translate -ot uint16 -srcwin $x $y $w $h data/PAN/img_$im.ntf tmp/data/crop/crop_img_$im.tif

qauto tmp/data/crop/crop_img_$im.tif tmp/data/crop/crop_img_$im.png

bin/create_mesh data/mcdsm/outdir_$im*/mcdsm.tif -21 tmp/output/small_mesh_$im.off tmp/output/small_mesh_$im.ply

bin/zbuffer \
    data/mcdsm/outdir_$im*/mcdsm.tif \
    -21 \
    tmp/data/crop/crop_img_$im.tif  \
    data/PAN/img_$im.rpc \
    tmp/output/xywh/xywh_$im.txt \
    tmp/output/small_mesh_$im
    tmp/output/matches/matches_img_$im.tif \

bin/colormultiple_mesh \
    tmp/output/small_mesh_$im.ply \
    tmp/data/crop/crop_img_$im.png \
    tmp/output/matches/matches_img_$im.tif  \
    tmp/output/coloured_$im.ply
