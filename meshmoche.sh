#!/binbash
# auteur : mariedautume

# exemple input bash meshmoche.sh /home/coco/src/s2p-iarpa/input/ essai_small/ data/small_lidar.tif
set -e
set -x
img_dir=$4
exp_dir=$1
lidar=$2
rpc_dir=$3
ox=$5
oy=$6

#IDX=`echo {01..42} {44..47}`
IDX=$7
mkdir -p $exp_dir/data/
mkdir -p $exp_dir/data/minMax/
#mkdir -p $exp_dir/data/cropped/

# crée à partir des grandes images .ntf des petites images tif contenant la zone du lidar 
echo "GET_CORNERS: get lidar corners projections on each image"
for i in $IDX; do
    bin/get_corners $lidar -21 $rpc_dir/img_$i.rpc > $exp_dir/data/minMax/minMaxWH_img_$i.txt
#    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
#    read ox oy < data/msi/offset/msi_offset_img_$i.txt
##   echo "bash colors/croprgb.sh $xmin $ymin $width $height $ox $oy $img_dir/img_$i.ntf data/msi/img_$i.ntf $exp_dir/data/cropped/cropped_img_$i.tif" 
#    bash colors/croprgb.sh $xmin $ymin $width $height $ox $oy $img_dir/img_$i.ntf data/msi/img_$i.ntf $exp_dir/data/cropped/cropped_img_$i.tif 
##    gdal_translate -ot uint16 -srcwin $xmin $ymin $width $height $img_dir/img_$i.ntf $exp_dir/data/cropped/cropped_img_$i.tif
done
#
#mkdir -p $exp_dir/data/images/

# quantifie les petites images tif et les convertit en jpg
#for i in $exp_dir/data/cropped/*tif; do
#    qeasy 100 1100 $i | plambda - "x 0.9 * x x join3" -o ${i/tif/png} 
#    convert ${i/tif/png} ${i/tif/jpg}
#done
#mv $exp_dir/data/cropped/*jpg $exp_dir/data/images/

echo $IDX

# construit la matrice de projection pour chaque image
mkdir -p $exp_dir/data/proj
echo "GET_P_OF_CROP: get projection matrix for each view"

for i in $IDX; do 
    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
    bin/get_P_of_crop $lidar -21 $rpc_dir/img_$i.rpc $xmin $ymin  > $exp_dir/data/proj/P_img_$i.txt
done

## calcule K, R et C à partir de P
#for i in $IDX; do 
#    python src/decomp_affine.py $exp_dir/data/proj/P_img_$i.txt
#done
#
mkdir -p $exp_dir/data/matches

# localise le lidar sur l'image puis reprojette sur lidar (donne les images avec grille)
echo "COLORIZE_WITH_SHADOWS: get lidar texture coordinates in each view"
for i in $IDX; do
    bin/zbuffer $lidar $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/proj/P_img_$i.txt $exp_dir/data/matches/matches_lidar_img_$i.tif
done
#
# crée un atlas et un mesh texturé à partir du lidar et d'une image.
echo "COLORSINGLE: create colorized mesh from lidar and one view"
mkdir -p $exp_dir/mesh
for i in $IDX; do
    bin/colorsingle $lidar $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/matches/matches_lidar_img_$i.tif $exp_dir/mesh/pil_$i.ply $exp_dir/mesh/atlas_$i -ox $ox -oy $oy
    qeasy 100 1100 $exp_dir/mesh/atlas_$i.tif $exp_dir/mesh/atlas_$i.png
done
##
## création d'un mesh à partir des 9 premières images
#echo "COLORMULTIPLE: create textured mesh from lidar and several views"
#bin/colormultiple $lidar $exp_dir/data/cropped/cropped_img_0*.tif data/rpc/img_0*.rpc $exp_dir/data/matches/matches_lidar_img_0*.tif $exp_dir/mesh/pil_multi.ply $exp_dir/mesh/atlas_multi 
#
#qeasy 100 1100 $exp_dir/mesh/atlas_multi.tif $exp_dir/mesh/atlas_multi.png
##
##  Created by marie d'autume on 10/05/2017.
##
