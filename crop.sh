#!/bin/bash
# auteur : mariedautume

img_dir=$1
exp_dir=$2

#mkdir $exp_dir/
#mkdir $exp_dir/data/
#mkdir $exp_dir/data/minMax/
#mkdir $exp_dir/data/cropped/
#
#for i in {01..47}; do
#    bin/get_corners data/Challenge1_Lidar.tif -21 $img_dir/img_$i.rpc > $exp_dir/data/minMax/minMaxWH_img_$i.txt
#    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
#    gdal_translate -ot uint16 -srcwin $xmin $ymin $width $height $img_dir/img_$i.ntf $exp_dir/data/cropped/cropped_img_$i.tif
#done

#mkdir $exp_dir/images/
#
#for i in $exp_dir/cropped/*tif; do
#    qeasy 100 1100 $i ${i/tif/png}
#    convert ${i/tif/png} ${i/tif/jpg}
#done
#mv $exp_dir/cropped/*jpg $exp_dir/cropped_jpg/
#

#mkdir $exp_dir/data/proj
#
#for i in {01..47}; do 
#    read xmin ymin width height < $exp_dir/data/minMax/minMaxWH_img_$i.txt
#    echo "bin/get_P_of_crop data/Challenge1_Lidar.tif -21 $exp_dir/data/rpc/img_$i.rpc $xmin $ymin essai.tif > $exp_dir/data/proj/P_img_$i.txt"
#    bin/get_P_of_crop data/Challenge1_Lidar.tif -21 $exp_dir/data/rpc/img_$i.rpc $xmin $ymin essai.tif > $exp_dir/data/proj/P_img_$i.txt
#done
#
#for i in {01..47}; do 
#    python src/decomp_affine.py essai_iarpa/data/proj/P_img_$i.txt
#done
mkdir $exp_dir/data/matches
for i in {01..47}; do
    echo "bin/colorize_with_shadows data/Challenge1_Lidar.tif $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/proj/P_img_$i.txt $exp_dir/data/matches/matches_lidar_img_$i.tif"
    bin/colorize_with_shadows data/Challenge1_Lidar.tif $exp_dir/data/cropped/cropped_img_$i.tif $exp_dir/data/proj/P_img_$i.txt $exp_dir/data/matches/matches_lidar_img_$i.tif
done
#
#
#  Created by marie d'autume on 10/05/2017.
#
