//#include "openMVG/sfm/sfm.hpp"
//#include "openMVG/image/image.hpp"

#include "Eigen"
#define _USE_EIGEN
#include <vector>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <string>

extern "C" {
#include "iio.h"
}
#include "../../openMVS/libs/MVS/Interface.h"
bool export_to_openMVS(int argc, char *argv[])
{
        std::string exp_dir = argv[1];
	char *filename_dsm = argv[2];

        MVS::Interface scene;
        size_t nPoses(0);
        const uint32_t nViews(5);
	int im_first = 0;

        for (int i=im_first; i<im_first+nViews; i++)
        {
                MVS::Interface::Platform platform;
                MVS::Interface::Platform::Camera camera;
                std::ifstream Kf;
                std::string Kfilename = exp_dir + "data/proj/K_img_0" + std::to_string(i+1) + ".txt";
                Kf.open(Kfilename);
                double x;
                for (int ci=0; ci<3; ci++)
                        for (int cj=0; cj<3; cj++)
                        {
                                Kf >> x;
                                std::cout << "impression de K :" << x << std::endl;
                                camera.K(ci,cj)=x;
                        }
                Kf.close();
                for (int c=0; c<3; c++)
                        camera.R(c,c) = 1.;
                platform.cameras.push_back(camera);
                scene.platforms.push_back(platform);
        }
        std::cout << scene.platforms[0].cameras.size() << std::endl;

        scene.images.reserve(nViews);
        for (int i=im_first; i<im_first+nViews; i++)
        {
                MVS::Interface::Image image;
                image.name = exp_dir + "data/images/cropped_img_0" + std::to_string(i+1) + ".jpg";
                image.platformID = i;
                MVS::Interface::Platform& platform = scene.platforms[image.platformID];
                image.cameraID = platform.cameras.size();
                MVS::Interface::Platform::Pose pose;
                image.poseID = platform.poses.size();
                std::ifstream Rf;
                std::string Rfilename = exp_dir + "data/proj/R_img_0" + std::to_string(i+1) + ".txt";
                Rf.open(Rfilename);
                double x;
                for (int ci=0; ci<3; ci++)
                        for (int cj=0; cj<3; cj++)
                        {
                                Rf >> x;
                                pose.R(ci,cj)=x;
                        }
                Rf.close();
                std::ifstream Cf;
                std::string Cfilename = exp_dir + "data/proj/R_img_0" + std::to_string(i+1) + ".txt";
                Cf.open(Cfilename);
                Cf >> x;
                pose.C.x = x;
                Cf >> x;
                pose.C.y = x;
                Cf >> x;
                pose.C.z = x;
                Cf.close();

                platform.poses.push_back(pose);
                std::cout << platform.poses.size() << std::endl;
                ++nPoses;
                scene.images.push_back(image);

        }
        std::cout << scene.platforms[0].poses.size() << std::endl;
        int w, h, pd;
        float *lidar = iio_read_image_float_vec(filename_dsm, &w, &h, &pd);

        int s = 10;
        scene.vertices.reserve(w*h/(s*s));
        for (int ci=0; ci<w/s; ci++)
        {
                std::cout << "nombre de colonnes sur total : " << 1.0*ci*s/w << std::endl;
                for (int cj=0; cj<h/s; cj++)
                        if (lidar[ci*s+w*cj*s]>0)
                        {
                                MVS::Interface::Vertex vert;
                                MVS::Interface::Vertex::ViewArr& views = vert.views;
        			for (int c=im_first; c<im_first+nViews; c++)
                                {
                                        std::string filename_match = exp_dir + "data/matches/matches_lidar_img_0" + std::to_string(c+1) + ".tif";
                                        const char * f_m = filename_match.c_str();
                                        int wm, hm, pdm;
                                        float *x = iio_read_image_float_vec(f_m, &wm, &hm, &pdm);
                                        if (!std::isnan(x[2*(cj*s*w+ci*s)]))
                                        {
                                                MVS::Interface::Vertex::View view;
                                                view.imageID = c;
                                                view.confidence = 0;
                                                views.push_back(view);
                                        }
                                        free(x);
                                }
                                vert.X.x = ci*s*1.0;
                                vert.X.y = cj*s*1.0;
                                vert.X.z = lidar[ci*s+w*cj*s];
                                scene.vertices.push_back(vert);
                        }
        }
        free(lidar);

	std::string filename_out = exp_dir + "scene.mvs";
        if (!MVS::ARCHIVE::SerializeSave(scene, filename_out))
                return false;

        std::cout
                << "Scene saved to openMVS" << std::endl;

        return true;
}

int main(int argc, char *argv[])
{
        export_to_openMVS(argc, argv);
        return 0;
}
