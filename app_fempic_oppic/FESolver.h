/*==============================================================================*
 * FESOLVER
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef FESOLVER_H
#define FESOLVER_H

#include <oppic_lib.h>
#include "fempic_ori/meshes.h"
#include <memory>

// #define USE_PETSC // define using makefile

#ifdef USE_PETSC
    #include <petscksp.h>
#endif

const double EPS0 = 8.8541878e-12;   /*permittivity of free space*/
const double QE   = 1.602e-19;       /*elementary charge*/
const double AMU  = 1.660538921e-27; /*atomic mass unit*/
const double Kb    = 8.617333262e-5;    /*Boltzmann's  constant*/

/*solver class*/
class FESolver {
public:
    enum Method {NonLinear, GaussSeidel, Lapack, Petsc};
    // double **K;        /*global stiffness matrix, should use a sparse matrix*/
    // double **J;        /*Jacobian matrix*/
    // double *F0;        /*"fh" and "fg" parts of the global force vector*/
    // double *F1;        /*"ff" part of the force vector*/

    int *ID;        /*ID[n]=A*/
    int **LM;        /*LM[e][a] location matrix */
    double ***NX;    /*NX[e][a] is a dNa/dx [3] vector*/
    int neq;        /*number of unknowns/equations*/

    /*reference values for the Boltzmann term*/
    double n0;
    double phi0;
    double kTe;
    double wall_potential;

    /*solution*/
    double *d;        /*d[neq] is the solution on the uknown nodes*/
    double *g;        /*g[n] essential boundaries*/
    double *uh;        /*uh[n] solution on nodes, union of d and g*/

    // double **ef;    /*ef[e][3] is the electric field in cell e*/

    double *detJ; /*determinant of the jacobian x_xi*/

    FESolver(opp::Params& params, std::shared_ptr<Volume> volume, int argc, char **argv);    /*constructor, initialized data structures*/
    ~FESolver();    /*destructor, frees memory*/

    void startAssembly();    /*clears K and F*/
    void preAssembly();
    void addKe(double** K, int e, double ke[4][4]);    /*adds contributions from element stiffness matrix*/
    // void addFe(double *F, int e, double fe[4]); /*adds contributions from element force vector*/
    void addFe(Vec *Fvec, int e, double fe[4]);

    double evalNa(int a, double xi, double eta, double zeta);
    void getNax(double nx[3], int e, int a);
    void inverse(double M[3][3], double V[3][3]);
    void computePhi(Method method, oppic_arg arg0, oppic_arg arg1);
    void buildF1Vector(double *ion_den);
    void solve(double *d, Method method);
    void solveNonLinear(double *d, double *y, double *G);
    void solveLinear(double **K, double *d, double *F);    /*solves Kd=F for d*/
    void solveLinearLapack(double **K, double *d, double *F);    /*solves Kd=F for d*/
    void initialzeMatrix(double **p_A);
    void solveLinearPetsc(double **K, double *d, double *F);  
    
    void updateEf();

    /*evaluates ef in cell e. Since constant field in cell, just copy*/
    // void evalEf(double res[3], int e) {for (int i=0;i<3;i++) res[i]=ef[e][i];}

    void summarize(std::ostream &out);

    void buildJmatrix(Method method);

    inline void setPotentialArray(double* potential_array) { uh = potential_array; }

protected:
    void computeNX();

    std::shared_ptr<Volume> volume;
    int n_nodes;
    int n_elements;    /*save this so we can properly deallocate LM*/

    /*quadrature points*/
    double l[2];
    double W[2];
    int n_int;

#ifdef USE_PETSC
    /* Petsc related variables */
    Vec         Xvec, Bvec, F0vec, F1vec, Dvec, Yvec;       
    Mat         Jmat, Kmat;                                 
    KSP         ksp;                                        /* linear solver context */
    KSPConvergedReason reason;

    int *vecCol, *matCol;
    int **matIndex;
    double **matIndexValues;

    bool matIndexCreated = false;
#endif
};


#endif /* !FESOLVER_H */
