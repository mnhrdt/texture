#:/usr/bin/env python

import sys
import os
import json
import numpy as np
import glob
import piio
from PIL import Image
import math
import tifffile

exp_dir = sys.argv[1]
print exp_dir
im_dir = os.path.join(exp_dir,'data/images/')
mat_dir = os.path.join(exp_dir,'data/proj')

im_path_list = glob.glob(os.path.join(im_dir,'*.jpg'))

sfm_data = {"sfm_data_version": "0.3",
        "root_path": im_dir,
        "views": [
            ],
        "intrinsics": [
            ],
        "extrinsics": [
            ],
        "structure": [
            ],
        "control_points": []
        };

#nb_images = len(im_path_list)
nb_images = 10

for i in np.arange(nb_images):

    im_path = im_path_list[i]
    im_name = os.path.split(im_path)[1]
    im = Image.open(im_path)
    im_number = im_name.split('.')[0].split('img_')[-1]
    
    sfm_data["views"].append({
        "key":i,
        "value":{
            "polymorphic_id": 0,
            "ptr_wrapper":{
                "data": {
                    "local_path": "",
                    "filename": im_name,
                    "width": im.size[0],
                    "height": im.size[1],
                    "id_view": i,
                    "id_intrisics": i,
                    "id_pose": i
                    }
                }
            }
        })

    K = np.loadtxt(os.path.join(mat_dir,'K_img_'+im_number+'.txt'))

    sfm_data["intrinsics"].append({
        "key":i,
        "value": {
            "polymorphic_id": 0,
            "ptr_wrapper": {
                "data": {
                    "width": im.size[0],
                    "height": im.size[1],
                    "focal_length": K[0,0],
                    "principal_point": [
                        K[0,2],
                        K[1,2]
                        ],
                    "disto_k3": [
                        0.0,
                        0.0,
                        0.0
                        ]
                    }
                }
            }
        })
    P = np.loadtxt(os.path.join(mat_dir,'P_img_'+im_number+'.txt'));
    R = np.loadtxt(os.path.join(mat_dir,'R_img_'+im_number+'.txt'));
    C = np.loadtxt(os.path.join(mat_dir,'C_img_'+im_number+'.txt'));

    sfm_data["extrinsics"].append({
        "key": i,
        "value": {
            "rotation": [
                [
                    R[0,0],
                    R[0,1],
                    R[0,2],
                    ],
                [
                    R[1,0],
                    R[1,1],
                    R[1,2],
                    ],
                [
                    R[2,0],
                    R[2,1],
                    R[2,2],
                    ]
                ],
            "center": [
                C[0],
                C[1],
                C[2]
                ]
            }
        })

lidar = tifffile.imread('data/Challenge1_Lidar.tif')
h = len(lidar)
w = len(lidar[0])

npoints = 0
s = 200
for i in np.arange(h/s):
    for j in np.arange(w/s):
        if lidar[i*s][j*s]>0.:
            sfm_data["structure"].append({
                "key": npoints,
                "value": {
                    "X": [
                        j*s,
                        i*s,
                        float(lidar[i*s][j*s])
                        ],
                    "observations":[]
                    }
                })
            for k in np.arange(nb_images):
                im_path = im_path_list[k]
                im_name = os.path.split(im_path)[1]
                im_number = im_name.split('.')[0].split('img_')[-1]

                match = tifffile.imread(os.path.join(exp_dir,'data/matches/matches_lidar_img_'+im_number+'.tif'))
           
                if not math.isnan(match[i*s][j*s][0]):
                    sfm_data["structure"][npoints]["value"]["observations"].append({
                        "key": k,
                        "value": {
                            "id_feat": 0,
                            "x": [
                                float(match[i*s][j*s][0]),
                                float(match[i*s][j*s][1])
                                ]
                            }
                        })
            npoints += 1
            print npoints
print npoints




with open(os.path.join(exp_dir,'data/sfm_data.json'),'w') as f:
    json.dump(sfm_data,f, indent=4)
#print json.dumps(sfm_data,indent=4)

                



    


    















