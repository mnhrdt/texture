#####
##### 3.5 Extract the internal parameters, rotation and center of the camera from the projection matrix
#####

# rq decomposition of M gives rotation and calibration
from scipy import linalg
import sys
import numpy as np
import os.path

P = np.loadtxt(sys.argv[1])
M = P[:,:-1]
M[-1,-1]=1.0
K, R = linalg.rq(M)
# fix sign of the scale params
R = np.diag(np.sign(np.diag(K))).dot(R)
K = K.dot(np.diag(np.sign(np.diag(K))))
R[-1,-1]=0.0

  # center
C = np.array([[P[0,1]*P[1,2]-P[1,1]*P[0,2]],
    [-P[0,0]*P[1,2]+P[1,0]*P[0,2]],
    [P[0,0]*P[1,1]-P[1,0]*P[0,1]],
    [0.0]])  

np.savetxt(sys.argv[1].replace('P','K'),K,fmt='%.20f')
np.savetxt(sys.argv[1].replace('P','R'),R,fmt='%.20f')
np.savetxt(sys.argv[1].replace('P','C'),C,fmt='%.20f')

