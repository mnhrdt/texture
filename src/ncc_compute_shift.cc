/* Copyright 2014, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr> */
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
extern "C" {
#include "iio.h"
}


//// a structure to wrap images
#include "img.h"
#include "img_tools.h"


//   struct Img corr(range*2+1,range*2+1);
void mean_var(struct Img &u, struct Img &msk, double *mu, double *var) {
   double muu=0, sigu=0;
   int count=0;
//#pragma omp parallel for reduction(+:count) reduction(+:muu)
   for (int y=0; y<u.ny; y++)
   for (int x=0; x<u.nx; x++) {
      float vu = valnan(u, x, y);
      if(valnan(msk,x,y)>0 && isfinite(vu)) {
         muu+=vu;
         count++;
      }
   }
   muu/=count;

   //var
//#pragma omp parallel for reduction(+:sigu)
   for (int y=0; y<u.ny; y++)
   for (int x=0; x<u.nx; x++) {
      float vu = valnan(u, x, y)-muu;
      if(valnan(msk,x,y) && isfinite(vu)) {
         sigu+=(vu*vu);
      }
   }
   sigu=sqrt(sigu/count);

   *var = sigu;
   *mu  = muu;
}


float ncc(struct Img &u, struct Img &v, int dx=0, int dy=0) {
   double muu=0, muv=0, sigu=0, sigv=0, xcorr=0;
   int count=0;
   //mean
//#pragma omp parallel for reduction(+:count) reduction(+:muu) reduction(+:muv)
   for (int y=0; y<u.ny; y++)
   for (int x=0; x<u.nx; x++) {
      float vu = valnan(u, x, y);
      float vv = valnan(v, x+dx, y+dy);
      if(isfinite(vu) && isfinite(vv)) {
         muu+=vu;
         muv+=vv;
         count++;
      }
   }
   muu/=count;
   muv/=count;

   //var
//#pragma omp parallel for reduction(+:xcorr) reduction(+:sigu) reduction(+:sigv)
   for (int y=0; y<u.ny; y++)
   for (int x=0; x<u.nx; x++) {
      float vu = valnan(u, x, y)-muu;
      float vv = valnan(v, x+dx, y+dy)-muv;
      if(isfinite(vu) && isfinite(vv)) {
         sigu+=(vu*vu);
         sigv+=(vv*vv);
         xcorr+=(vu*vv);
      }
   }
   sigu=sqrt(sigu/count);
   sigv=sqrt(sigv/count);
   xcorr=xcorr/count;

   return  xcorr/(sigu*sigv);
}

void compute_robust_statistics(struct Img &u, struct Img &msk, double *mu, double *sig) {
   *mu=0; *sig=1;
   std::vector<float> tmp;
   for(int i=0;i<u.nx*u.ny;i++) if( isfinite(msk[i]) && isfinite(u[i]) ) tmp.push_back(u[i]);
   if(tmp.size()==0) return;
   std::sort(tmp.begin(), tmp.end()); //, std::greater<float>());
   *mu = tmp[tmp.size()/2];
//   *sig = tmp[(tmp.size()*3)/4] - tmp[tmp.size()/4];  // IQR
   for(int i=0;i<tmp.size();i++) tmp[i] = fabs(tmp[i]-*mu);
   std::sort(tmp.begin(), tmp.end()); //, std::greater<float>());
   *sig = 1.4826 * tmp[tmp.size()/2];                   //MAD
}


void compute_ncc(struct Img &u, struct Img &v, int range, int *dx, int *dy) {
   int initdx=*dx, initdy=*dy;
   float maxv=-INFINITY;
#pragma omp parallel for
   for (int y=initdy-range; y<=initdy+range; y++)
   for (int x=initdx-range; x<=initdx+range; x++) {
      float corr;
      corr = ncc(u, v, x, y);
#pragma omp critical
      if (corr > maxv) {
         *dx=x;
         *dy=y;
         maxv=corr;
      }
   }
}

// downsampling 2x of the image u
inline struct Img downsample2x(struct Img &u) {
   // allocate output
   struct Img o(ceil(u.nx/2.0),ceil(u.ny/2.0),u.nch);

   // apply the filter
   for(int c=0;c<u.nch;c++)
   for(int j=0;j<u.ny;j+=2)
   for(int i=0;i<u.nx;i+=2) {
      float v=0, count=0;
      for(int k=0;k<2;k++) for(int l=0;l<2;l++) {
         float t=valnan(u, Point(i+k,j+l));
         if (isfinite(t)) {
            v+=t;
            count++;
         }
      }
      o[i/2+j/2*o.nx + o.npix*c] = count>0 ? v/count: NAN;
   }
   return o;
}


void recursive_ncc(struct Img &u, struct Img &v, int range, int *dx, int *dy) {
   if(fmin(u.nx,u.ny) > 100) {
      struct Img su = downsample2x(u);
      struct Img sv = downsample2x(v);
      *dx /= 2; *dy /= 2;
      recursive_ncc(su, sv, range, dx, dy);
      *dx *= 2; *dy *= 2;
   }

   compute_ncc(u, v, range, dx, dy);
//   fflush(stdout);
//   printf("%d %d\n", *dx, *dy);
}




#ifndef SKIP_MAIN
// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
static char *pick_option(int *c, char ***v, char *o, char *d) {
   int argc = *c;
   char **argv = *v;
   int id = d ? 1 : 0;
   for (int i = 0; i < argc - id; i++)
      if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, o)) {
     char *r = argv[i + id] + 1 - id;
     *c -= id + 1;
     for (int j = i; j < argc - id; j++)
        (*v)[j] = (*v)[j + id + 1];
     return r;
      }
   return d;
}

int main(int argc, char* argv[])
{
    bool scaling = pick_option(&argc, &argv, "a", NULL);
    if (argc < 3) {
        fprintf(stderr, "too few parameters\n");
        fprintf(stderr, "\tusage: %s u v [range(5) [dx dy]]\n", argv[0]);
        return 1;
    }

    //read the parameters
    int i = 1;
    char* f_u     = (argc>i)? argv[i]      : (char*) "" ;      i++;
    char* f_v     = (argc>i)? argv[i]      : (char*) "" ;      i++;
    int   range   = atoi((argc>i)? argv[i] : (char*) "5");     i++;
    int   dx      = argc > i ? atoi(argv[i]) : 0;     i++;
    int   dy      = argc > i ? atoi(argv[i]) : 0;     i++;

    // read input
    struct Img u = iio_read_vector_split(f_u);
    struct Img v = iio_read_vector_split(f_v);
    struct Img outv(v);

    if (argc < 6)
        recursive_ncc(u, v, range, &dx, &dy);

     //iio_write_vector_split("corr.tif", corr);

    struct Img mskv(v.nx, v.ny);
    for (int i = 0; i < mskv.npix; i++)
        mskv[i] = 0;
    for (int y = 0; y < v.ny; y++)
    for (int x = 0; x < v.nx; x++) {
        outv(x, y) = valnan(v, x+dx, y+dy);
        float vu = valnan(u, x, y);
        float vv = valnan(outv, x, y);
        mskv(x, y) = (isfinite(vu) && isfinite(vv));
    }

    struct Img msk(u.nx, u.ny);
    for (int i = 0; i < msk.npix; i++)
        msk[i] = 0;
    for (int y = 0; y < u.ny; y++)
    for (int x = 0; x < u.nx; x++) {
        float vu = valnan(u, x, y);
        float vv = valnan(outv, x, y);
        msk(x,y) = (isfinite(vu) && isfinite(vv));
    }

    double muu, sigu, muv, sigv;
    mean_var(u,    msk , &muu, &sigu);
    mean_var(outv, mskv, &muv, &sigv);

    // delta x, delta y, delta v (affine as: alpha v + b)
    double a = scaling ? sigu/sigv : 1;
    double b = muu - muv*a;
    printf("%d %d %f %f\n", dx, dy, a, b);
    return 0;
}
#endif
