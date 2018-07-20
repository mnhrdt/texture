SCALE_X=0.300000011920929 # metres
SCALE_Y=-0.300000011920929 # metres
ORIG_X=354052.375000000000000 # easting-northing (metres)
ORIG_Y=6182702.000000000000000 # easting-northing (metres)
ZONE=-21 # signed utm zone
CROP_X=1505 # pixels
CROP_Y=891 # pixels
CROP_WIDTH=390 # pixels
CROP_HEIGHT=372
RES_MESH=0.3 # metres
#RECALAGE=bin/ncc_compute_shift 
RECALAGE=bin/gc
DSM=data/Challenge1_Lidar_nan.tif
INPUTS=data/input_try.txt

# ZONE=$(if [ $(gdalinfo $DSM | grep UTM | awk -F"zone" '{print $2}' | cut -f 1 -d \" | tail -c 2) = "S" ]; then echo $(gdalinfo $DSM | grep UTM | awk -F"zone" '{print $2}' | awk -F$(gdalinfo $DSM | grep UTM | awk -F"zone" '{print $2}' | cut -f 1 -d \" | tail -c 2) '{print -$1}'); else echo $(gdalinfo $DSM | grep UTM | awk -F"zone" '{print $2}' | awk -F$(gdalinfo $DSM | grep UTM | awk -F"zone" '{print $2}' | cut -f 1 -d \" | tail -c 2)  '{print $1}'); fi);

# SCALE_X=$(gdalinfo $DSM | grep Pixel | cut -f 2 -d \( | cut -f 1 -d ,
# SCALE_X=$(gdalinfo $DSM | grep Pixel | cut -f 1 -d \) | cut -f 2 -d ,


# ORIG_X=$(gdalinfo $DSM | grep Origin | cut -f 2 -d \( | cut -f 1 -d ,
# ORIG_X=$(gdalinfo $DSM | grep Origin | cut -f 1 -d \) | cut -f 2 -d ,

# RES_MESH=$SCALE_X


MESH=output_try/mesh/refined_mesh.off
