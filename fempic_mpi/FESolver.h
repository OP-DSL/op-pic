/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include <oppic_lib.h>
#include "fempic_defs.h"
#include "fempic_ori/meshes.h"
#include "fempic_ori/maths.h"

#ifndef USE_PETSC
    #define Vec int
    #define Mat int
    #define KSP int
    #define KSPConvergedReason int
#endif

const double EPS0 = 8.8541878e-12;      /*permittivity of free space*/
const double QE   = 1.602e-19;          /*elementary charge*/
const double AMU  = 1.660538921e-27;    /*atomic mass unit*/
const double Kb   = 8.617333262e-5;     /*Boltzmann's  constant*/

/*solver class*/
class FESolver {
public:
    enum Method { NonLinear, GaussSeidel, Lapack, Petsc };

    FESolver(
        oppic_map cell_to_nodes_map, 
        oppic_dat node_type, 
        oppic_dat node_pos,  
        oppic_dat node_bnd_pot,
        int argc, char **argv);
    ~FESolver();

    void computePhi(oppic_arg arg0, oppic_arg arg1, oppic_arg arg2);
    
    void preAssembly(oppic_map cell_to_nodes_map, oppic_dat node_bnd_pot);
    void enrich_cell_shape_deriv(oppic_dat cell_shape_deriv);

    bool linearSolve(double *ion_den); 
    void nonLinearSolve(double *ion_den); 
    void buildJmatrix();
    void buildF1Vector(double *ion_den);
    
    void summarize(std::ostream &out);  

protected:
    void addKe(double** K, int e, double ke[4][4]);
    void addFe(Vec *Fvec, int e, double fe[4]);
    double evalNa(int a, double xi, double eta, double zeta);
    void getNax(double nx[3], int e, int a);
    void initialzeMatrix(double **p_A);
    void computeNX(oppic_dat node_pos, oppic_map cell_to_nodes_map);
    void sanityCheck();
    void initID(oppic_dat node_type_dat);
    void initPetscStructures();

    Method fesolver_method = Method::Petsc;

    int *ID;        /*ID[n]=A*/
    int **LM;       /*LM[e][a] location matrix */
    double ***NX;   /*NX[e][a] is a dNa/dx [3] vector*/
    int neq;        /*number of unknowns/equations*/
    double *detJ;     /* determinant of the jacobian x_xi */

    double n0;
    double phi0;
    double kTe;
    double wall_potential;

    double *d;        /* d[neq] is the solution from the linear solver */  

    int n_nodes_set = 0;
    int n_nodes_inc_halo = 0; 
    int n_elements_set = 0;
    int n_elements_inc_halo = 0;

    int own_start = 0;
    int own_end = 0;
    int global_neq = 0;

    /*quadrature points*/
    double l[2];
    double W[2];
    int n_int = 0;

    /* Petsc related variables */
    Vec         Xvec, Bvec, F0vec, F1vec, Dvec, Yvec;       
    Mat         Jmat, Kmat;                                 
    KSP         ksp;            /* linear solver context */
    KSPConvergedReason reason;

    // TODO : revise -- these can be removed...
    int *vecCol;                    
    int *matCol;                    // Number of non zero columns per row
    int **matIndex;                 // Non zero column indices per row
    double **matIndexValues;        // Non zero column values per row (temp copy to enrich K matrix)
    bool matIndexCreated = false;   // Not required
};

