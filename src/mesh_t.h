#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>
#include <libgen.h>
#include <math.h>
#include <string.h>

#include "rpc.c"
#include "fail.c"
#include "iio.h"


// Création d'une structure pour le mesh. À optimiser

struct image_coord{ // coordonnées de la projection du sommet sur l'image
        double i; 
        double j;
};

struct vertex{
        int ij[2]; // coordinates in lidar
        double xyz[3]; // coord dans l'espace
        struct image_coord *im; // liste coordonnées sommet dans chaque im 
};

struct face{
        int v0; // 1er sommet
        int v1; // 2ème sommet
        int v2; // 3ème sommet
        double n[3]; // vecteur normal à la face
        int im; // référence de l'image la plus en face
        bool *isvisible; // 1 si visible dans l'image i, 0 sinon
};

struct mesh_t{
        int nimages; // nombre d'images
        int nv; // nombre de sommets
        int nf; // nombre de faces
        struct vertex *v; // liste des sommets
        struct face *f; // liste des faces
};



