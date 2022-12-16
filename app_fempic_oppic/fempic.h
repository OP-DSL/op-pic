// Using this file since the original fempic mesh loaders and functions are used

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>
#include <math.h>
#include <array>
#include <chrono>
#include <algorithm> 
#include "../lib_oppic/trace.h"

#define USE_PETSC
#define USE_PARTICLE_SORTING

#ifdef USE_PETSC
    #include <petscksp.h>
#endif

using namespace std;

#define DIMENSIONS         3
#define NODES_PER_CELL     4
#define NEIGHBOUR_CELLS    4
#define DET_FIELDS         4
#define PRINT_PRECISION    15

/*constants*/
const bool print_all_to_file = false;

const double EPS0     = 8.8541878e-12;    /*permittivity of free space*/
const double QE       = 1.602e-19;        /*cellary charge*/
const double AMU      = 1.660538921e-27;  // kg, atomic mass unit

const double PLASMA_DEN      = 1e10; 
const double ION_VELOCITY    = 7000;    

const double OP_CONST_charge = 1*QE;     // The ratio of real particles per macro-particle is called the specific weight,
const double OP_CONST_mass   = 16*AMU;
const double OP_CONST_spwt   = 2e2;

const bool bool_true         = true;
const bool bool_false        = false;

const double area            = 0.2*0.2;                            /*set area of the k=0 face, this should be the sum of triangle areas on the inlet*/       
const double num_per_sec     = PLASMA_DEN*ION_VELOCITY*area;       /*number of real ions per sec, given prescribed density and velocity*/     

/*node type*/
enum NodeType {NORMAL,OPEN,INLET,SPHERE};

/*definition of a node*/
struct Node
{
    Node(double x, double y, double z) 
    {
        pos[0]=x;
        pos[1]=y;
        pos[2]=z;
        type=NORMAL;
    }    
    double pos[3];    /*node position*/
    NodeType type;
    double volume;    /*node volume*/
};

/*definition of a tetrahedron*/
struct Tetra
{
    int con[NODES_PER_CELL];
    double volume;
    Tetra (int n1, int n2, int n3, int n4) 
    {
        con[0]=n1;
        con[1]=n2;
        con[2]=n3;
        con[3]=n4;
    }
    
    /*data structures to hold precomputed 3x3 determinants*/
    double alpha[4], beta[4], gamma[4], delta[4];
    
    /*cell connectivity*/
    int cell_con[NEIGHBOUR_CELLS];    /*index corresponds to the face opposite the i-th node*/    
};

/*definition of a volume*/
struct Volume
{
    vector <Node> nodes;
    vector <Tetra> cells;
};


/*species class*/
class Particles
{
public:
    vector<int> particle_cell_index;
    vector<array<double,DIMENSIONS>> particle_pos;
    vector<array<double,DIMENSIONS>> particle_vel;
    vector<array<double,DIMENSIONS>> particle_ef;
    vector<array<double,NODES_PER_CELL>> particle_lc;    /*particle's weights*/

    vector<int> particle_indexes_to_remove;
    vector<int> particle_cell_index_to_add;
    vector<array<double,DIMENSIONS>> particle_pos_to_add;
    vector<array<double,DIMENSIONS>> particle_vel_to_add;
    vector<array<double,DIMENSIONS>> particle_ef_to_add;
    vector<array<double,NODES_PER_CELL>> particle_lc_to_add;    /*particle's weights*/

    double spwt; // The ratio of real particles per macro-particle is called the specific weight,
    double mass;
    double charge;
    double rem;        /*fractional particle reminder*/

    int get_num_particles();

    void collect_particle_to_add(
        const int cell_index, 
        array<double,DIMENSIONS> pos, 
        array<double,DIMENSIONS> vel, 
        array<double,DIMENSIONS> ef, 
        array<double,NODES_PER_CELL> lc);

    void insert_collected_particles();

    void mark_to_remove_particle(const int particle_index);

    void remove_marked_particles();

    Particles () 
    {
        rem=0;
    }
    ~Particles () 
    {

    }
};

/*solver class*/
class FESolver
{
public:
    double **p_K;      /*global stiffness matrix, should use a sparse matrix*/
    double **p_J;      /*Jacobian matrix*/
    double *F0;        /*"fh" and "fg" parts of the global force vector*/
    double *F1;        /*"ff" part of the force vector*/

    int *ID;        /*ID[n]=A*/
    int **LM;       /*LM[e][a] location matrix */
    double ***NX;   /*NX[e][a] is a dNa/dx [3] vector*/
    int neq;        /*number of unknowns/equations*/

    /*reference values for the Boltzmann term*/
    double n0;
    double phi0;
    double kTe;

    /*solution*/
    double *p_d;    /*d[neq] is the solution on the uknown nodes*/

    double *detJ; /*determinant of the jacobian x_xi*/    
    
    FESolver(Volume &volume);    /*constructor, initialized data structures*/
    ~FESolver();    /*destructor, frees memory*/

    void SolveFields(
        const double *ion_den, 
        double *field_potential, 
        double *boundary_potential, 
        double *electric_field);

    void startAssembly();    /*clears K and F*/
    void preAssembly(double *boundary_potential);

private:
    void addKe(int e, double ke[4][4]);    /*adds contributions from cell stiffness matrix*/
    void addFe(double *F, int e, double fe[4]); /*adds contributions from cell force vector*/

    double evalNa(int a, double xi, double eta, double zeta);
    void getNax(double nx[3], int e, int a, double xi, double eta, double zeta);
    void inverse(double M[3][3], double V[3][3]);
    void computePhi(const double *ion_den, double *field_potential, double *boundary_potential);    
    void buildF1Vector(const double *ion_den);
    void solveNonLinear(const double *ion_den);
    void solveLinear(double **K, double *d, double *F);    /*solves Kd=F for d*/
#ifdef USE_PETSC
    void initialze_matrix(double **p_A);
    void solveLinear_petsc(double **K, double *d, double *F);
#endif
    void updateElementElectricField(double *field_potential, double *electric_field);
    void buildJmatrix();
    void computeNX();
    
    Volume &volume;    
    int n_nodes;
    int n_cells;    /*save this so we can properly deallocate LM*/

    /*quadrature points*/
    double p_l[2];
    double p_W[2];
    int n_int;

#ifdef USE_PETSC
    /* Petsc related variables */
    Vec         x, b;        /* approx solution, RHS, exact solution */
    Mat         A;              /* linear system matrix */
    KSP         ksp;            /* linear solver context */

    int *vecCol, *matCol;
    int **matIndex;
    double **matIndexValues;

    bool matIndexCreated = false;
#endif
};


bool LoadVolumeMesh(const string file_name, Volume &volume);
bool LoadSurfaceMesh(const string file_name, Volume &volume, NodeType node_type);
void OutputMesh(int ts, Volume &volume, double *phi, double *ef, double *ion_den);
void OutputParticles(int ts, Particles &particles);
double det4(double (*M)[4]);
double det3(double (*M)[3]);
void matVecMultiply(double *y, double**A, double *x, int nu);
void vecVecSubtract(double *, double *v1, double *v2, int nu);
bool XtoLtet(int& part_cell_index, array<double, 4>& part_lc, const array<double, DIMENSIONS>& part_pos, const Volume &volume, bool search=true);
double rnd();

void print_particles_m(const double* part_pos, const double* part_vel, const int* cell_index, int n_particles, const std::string suffix);
void print_fields_m(double *ion_den, double *field_potential, double *electric_field, int n_nodes, int n_cells, const std::string suffix);
void EnrichArrays(
    Volume& volume, 
    double *&node_bnd_pot_tmp, double *&node_pot_tmp, double *&node_pos_tmp, double *&node_volume_tmp, double *&node_ion_den_tmp,
    int *&cell_to_nodes_tmp, int *&cell_to_cell_tmp, double *&cell_det_tmp, double *&cell_volume_tmp, double *&cell_ef_tmp);
int get_num_particles_to_inject(double dt, double& remainder);
