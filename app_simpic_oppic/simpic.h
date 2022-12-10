#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <stdarg.h>

// global variables
#define MPI_double MPI_DOUBLE
#define STR_SIZE 		501
#define NG 				100
#define NP 				10

#define DEBUG

// definitions for diagnostics
#define NDENSITY		0
#define NEFIELD 		1
#define NVOLTAGE		2
#define NVXPHASE		3
#define NDIAGNOSTICS 	4

// variables for timing
#define PARTICLES		0
#define FIELDS			1
#define MPITIME			2
#define DIAGNOSTICS		3
#define INIT			4
#define NTIMES			5

const double fourth			 	= 1.0 / 4.0;
const double half			   	= 1.0 / 2.0;
const double one				= 1.;
const double one_third		  	= 1./3.;
const double two_fifteenths	 	= 2./15.;

extern double dtfactor;

extern double *pdata;
extern double *lhsbuf, *rhsbuf;
extern int npart, lpart, rpart;

// fields arrays
extern double *Earray;
extern double *phiarray;
extern double *narray;
extern double *nback;

extern double area;
extern double density;
extern double np2c;
extern double q, m, qm;
extern double qscale;
extern double epsilon;
extern double wp;
extern double El, Er;
extern double dx, dt, t;
extern double xl, xr, L, Ll;
extern int ntimesteps;

extern int lproc, rproc;
extern int nproc, rank;
extern int ver, subver;
extern int last;

extern int ng, nc, ngglobal;

extern int nl, nr;
extern int nlp, nrp;  // counters for sending particle buffers
extern int nlr, nrr;

extern double *tri_a, *tri_b, *tri_c, *gam; // TODO : Enrich properly

extern bool diagnosticsflag;
extern double lhsvoltage;

extern double tttimes[NTIMES];
extern char names[NTIMES][STR_SIZE];
extern double tstart, tend;  // global times
extern double wtres;
extern double nparttot;

extern double tt;  // universal time

// variables for field communication
extern double *frbuffer;
extern double *fsbuffer;

extern char diags[NDIAGNOSTICS][STR_SIZE];
extern double * thediags[NDIAGNOSTICS];

extern char procname[STR_SIZE];

//*************************************************************************************************
// global functions
void Quit(void);

inline double xcoord(int i);

void gradient(double *grad, double *arr, int n, double scale);

void parsecmdline(int argc, char **argv);

void init(void);

double rhsV(double t);

double lhsV(double t);

void setcoeffs(double scale);


//*************************************************************************************************
void init_particles(
	double *&particle_position_x_tmp, 
	double *&particle_velocity_x_tmp, 
	double *&particle_field_E_tmp, 
	int *&particle_cell_index_tmp);
	
void init_fields(
	double *&node_field_E_tmp, 
	double *&node_field_J_tmp, 
	double *&node_field_P_tmp,
	double *&node_xlocal_tmp, 
	int *&node_index_tmp, 
	int *&cell_to_nodes_tmp);

void printDataToFiles(
	int n_nodes,
	int n_particles,
	double* node_field_E, 
	double* node_field_J, 
	double* node_field_P,
	double* particle_position_x, 
	double* particle_velocity_x, 
	int* particle_cell_index,
	const std::string suffix);

//*************************************************************************************************
void seq_field_solve_poissons_equation(
	int set_size,
	double* node0_field_J,
	double* node0_field_P);

//*************************************************************************************************
void seq_field_solve_get_potential_gradient(
	int set_size,
	double* nodes_field_E,
	const double* nodes_field_P);

//*************************************************************************************************


