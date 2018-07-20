#!/bin/bash
# auteur : mariedautume

# exemple d'utilisation à partir de ../
# bash entire.sh config.sh contrast_method fusion_method output_dir
# exemple de config.sh dans scripts/config_try.sh
# contrast_method = mean_var pour l'instant
# fusion_method   = frontal_{vertices, edges}, scalar_{vertices, edges}

set -e
set -x
config=$1
c=$2
f=$3
output=$4

. $config

dir=`mktemp -d`

# create all sub folders
mkdir -p $dir/mesh
mkdir -p $dir/fused
mkdir -p $dir/ncc_shift
mkdir -p $output/mesh

# if no mesh provided in config.sh, create mesh from reference dsm
if [ -z "$MESH" ]
then

    gdal_translate -ot float64 -srcwin 1505 891 390 372 $DSM $dir/mesh/small_dsm.tif

    bin/create_mesh \
        $dir/mesh/small_dsm.tif \
        $dir/mesh/scaled_mesh.off \
        $dir/mesh/scaled_mesh.ply 

    bin/refine \
        $dir/mesh/scaled_mesh.off \
        $output/mesh/refined_mesh.off \
        --res $RES_MESH

    echo MESH=$output/mesh/refined_mesh.off >> $config
    MESH=$output/mesh/refined_mesh.off
fi

############ TO DO #################
# if mesh provided but no reference dsm create it from mesh
###################################

# get edges list from mesh
bin/triproc off2edges \
    $MESH \
    $dir/edges.txt

# loop over all images to get shadow and colour information for each vertex 
# from one image
cat $INPUTS | xargs -L1 ./scripts/one_image.sh $dir $config


############# FUSION #############
# all matlab files are supposed to be in a subfolder called scripts
folder=`pwd`

$OCTAVE "cd '$folder/scripts/'; matlab_wrapper('$dir', '$c', '$f', '$dir/edges.txt', '$dir/fused/fused_$c$f.tif'); exit;"

######## QUANTIFICATION #########

######## TO DO #################
# trouver une meilleure méthode
################################

qauto $dir/fused/fused_$c$f.tif $dir/fused/fused_$c\_$f.tiff

# création du mesh final
bin/write_coloured_ply \
    $MESH \
    $dir/fused/fused_$c\_$f.tiff \
    $output/mesh/fused_mesh_$c\_$f.ply
     

# suppression des dossiers temporaires
#rm -r $dir
