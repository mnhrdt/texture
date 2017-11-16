im_dir=$1
exp_dir=$2
focal=$3 # exemple 1663875 44.37

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_SfMInit_ImageListing -i $im_dir -f $focal -o $exp_dir/matches/

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_ComputeFeatures -i $exp_dir/matches/sfm_data.json -o $exp_dir/matches/ -m SIFT -f 1 -p ULTRA

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_ComputeMatches -i $exp_dir/matches/sfm_data.json -o $exp_dir/matches/ -f 1 -r .8  -g e

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_GlobalSfM -i $exp_dir/matches/sfm_data.json -m $exp_dir/matches/ -o $exp_dir/SfM

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_ComputeStructureFromKnownPoses  -i $exp_dir/SfM/sfm_data.bin -m $exp_dir/matches/ -o $exp_dir/SfM/robust.bin

/home/mariedautume/src/openMVG_build/Linux-x86_64-RELEASE/openMVG_main_openMVG2openMVS -i $exp_dir/SfM/sfm_data.bin -o $exp_dir/mvs/scene.mvs -d $exp_dir/mvs

/home/mariedautume/src/openMVS_build/bin/DensifyPointCloud $exp_dir/mvs/scene.mvs

/home/mariedautume/src/openMVS_build/bin/ReconstructMesh $exp_dir/mvs/scene_dense.mvs

/home/mariedautume/src/openMVS_build/bin/RefineMesh $exp_dir/mvs/scene_dense_mesh.mvs

/home/mariedautume/src/openMVS_build/bin/TextureMesh $exp_dir/mvs/scene_dense_mesh.mvs

/home/mariedautume/src/openMVS_build/bin/TextureMesh $exp_dir/mvs/scene_dense_mesh_refine.mvs
















