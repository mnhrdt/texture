#!/binbash
# auteur : mariedautume

# exemple input bash all_new.sh 01
set -e
set -x
config=$1
c=$2
f=$3
output=$4

. $config

dir=`mktemp -d`

mkdir -p $dir/mesh
mkdir -p $dir/fused
mkdir -p $dir/ncc_shift
mkdir -p $output/mesh

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

bin/triproc off2edges \
    $MESH \
    $dir/edges.txt

echo $INPUTS

cat $INPUTS

cat $INPUTS | xargs -L1 ./scripts/one_image.sh $dir $config


############# FUSION #############
folder=`pwd`

matlab -nodesktop -nojvm -nosplash -r\
    "cd '$folder/scripts/'; matlab_wrapper('$dir', '$c', '$f', '$dir/edges.txt', '$dir/fused/fused_$c$f.tif'); exit;"

######## QUANTIFICATION #########

qauto $dir/fused/fused_$c$f.tif $dir/fused/fused_$c\_$f.tiff

bin/write_coloured_ply \
    $MESH \
    $dir/fused/fused_$c\_$f.tiff \
    $output/mesh/fused_mesh_$c\_$f.ply
     
#rm -r $dir
