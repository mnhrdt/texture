#!/bin/bash
# auteur : mariedautume

# exemple input bash crop.sh /home/coco/src/s2p-iarpa/input/ essai/ data/small_Challenge1_Lidar.tif
set -e
img_dir=$1
exp_dir=$2
lidar=$3

IDX=`echo {01..05}`

mkdir -p $exp_dir/
mkdir -p $exp_dir/data/
mkdir -p $exp_dir/data/minMax/
mkdir -p $exp_dir/data/cropped/

# crée à partir des grandes images .ntf des petites images tif contenant la zone du lidar 
for i in $IDX; do
    bin/get_corners $lidar -21 $img_dir/img_$i.rpc > $exp_dir/data/minMax/minMaxWH_img_$i.txt
    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
    gdal_translate -ot uint16 -srcwin $xmin $ymin $width $height $img_dir/img_$i.ntf $exp_dir/data/cropped/cropped_img_$i.tif
done

mkdir -p $exp_dir/data/images/

# quantifie les petites images tif et les convertit en jpg
for i in $exp_dir/data/cropped/*tif; do
    qeasy 100 1100 $i | plambda - "x 0.9 * x x join3" -o ${i/tif/png} 
    convert ${i/tif/png} ${i/tif/jpg}
done
mv $exp_dir/data/cropped/*jpg $exp_dir/data/images/



mkdir -p $exp_dir/data/proj

# construit la matrice de projection pour chaque image
for i in $IDX; do 
    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
    bin/get_P_of_crop $lidar -21 $img_dir/img_$i.rpc $xmin $ymin essai.tif > $exp_dir/data/proj/P_img_$i.txt
done

# calcule K, R et C à partir de P
for i in $IDX; do 
    python src/decomp_affine.py $exp_dir/data/proj/P_img_$i.txt
done

mkdir -p $exp_dir/data/matches

# localise le lidar sur l'image puis reprojette sur lidar (donne les images avec grille)
for i in $IDX; do
    echo "bin/colorize_with_shadows data/Challenge1_Lidar.tif $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/proj/P_img_$i.txt $exp_dir/data/matches/matches_lidar_img_$i.tif"
    bin/colorize_with_shadows data/Challenge1_Lidar.tif $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/proj/P_img_$i.txt $exp_dir/data/matches/matches_lidar_img_$i.tif
done
#
#
#  Created by marie d'autume on 10/05/2017.
#
