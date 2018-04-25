#ifndef RPC_H
#define RPC_H
// rational polynomial coefficient stuff

#include <stdio.h>



// rational polynomial coefficients (and those of the inverse model)
struct rpc 
{
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];

	double dmval[4];
	double imval[4];
};

void read_rpc_file_ikonos(struct rpc *p, char *filename);

void read_rpc_file_xml_pleiades(struct rpc *p, char *filename);

void read_rpc_file_xml_worldview(struct rpc *p, char *filename);

// read a file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename);

#define FORI(n) for (int i = 0; i < (n); i++)

void print_rpc(FILE *f, struct rpc *p, char *n);

// evaluate a polynomial of degree 3
double eval_pol20(double c[20], double x, double y, double z);

double eval_pol20_dx(double c[20], double x, double y, double z);

double eval_pol20_dy(double c[20], double x, double y, double z);

double eval_pol20_dz(double c[20], double x, double y, double z);

// evaluate the direct rpc model
void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate the inverse rpc model
void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate a correspondence between two images given their rpc
void eval_rpc_pair(double xprime[2],
		struct rpc *pa, struct rpc *pb,
		double x, double y, double z);

void rpc_projection(double ij[2], struct rpc *r, double lonlatheight[3]);

void rpc_localization(double lonlat[2], struct rpc *r, double ijh[3]);


#define RPCH_MAXIT 100
#define RPCH_HSTEP 1
#define RPCH_LAMBDA_STOP 0.00001
#define RPCH_A2MAX 1e-50
// compute the height of a point given its location inside two images
double rpc_height(struct rpc *rpca, struct rpc *rpcb,
		double xa, double ya, double xb, double yb, double *outerr);


#endif
