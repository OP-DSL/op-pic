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

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "FESolver.h"
#include "fempic.h"

#define is_neq(a, b) ((a) != (b))

#define SCALLING 100

bool print_petsc = false;

//*************************************************************************************************
FESolver::FESolver(
    oppic_map cell_to_nodes_map, 
    oppic_dat node_type_dat, 
    oppic_dat node_pos_dat,
    oppic_dat node_bnd_pot_dat,
    int argc, char **argv) 
{
    if (OP_DEBUG) opp_printf("FESolver", "FESolver");

    phi0           = 0;
    n0             = opp_params->get<OPP_REAL>("plasma_den");
    kTe            = Kb * opp_params->get<OPP_REAL>("electron_temperature");
    wall_potential = -(opp_params->get<OPP_REAL>("wall_potential"));
    neq            = 0; /*count number of unknowns*/

    if      (opp_params->get<OPP_STRING>("fesolver_method") == "nonlinear") fesolver_method = NonLinear;
    else if (opp_params->get<OPP_STRING>("fesolver_method") == "gaussseidel") fesolver_method = GaussSeidel;
    else if (opp_params->get<OPP_STRING>("fesolver_method") == "lapack") fesolver_method = Lapack;
    else if (opp_params->get<OPP_STRING>("fesolver_method") == "petsc") fesolver_method = Petsc;  
    
    n_nodes_set         = node_type_dat->set->size;
    n_nodes_inc_halo    = node_type_dat->set->size + node_type_dat->set->exec_size + node_type_dat->set->nonexec_size;
    n_elements_set      = cell_to_nodes_map->from->size;
    n_elements_inc_halo = cell_to_nodes_map->from->size + cell_to_nodes_map->from->exec_size; 

    // opp_printf("FESolver", "n_nodes_set %d | n_nodes_inc_halo %d | n_elements_set %d | n_elements_inc_halo %d", 
    //     n_nodes_set, n_nodes_inc_halo, n_elements_set, n_elements_inc_halo);

    // TODO : DO NOT CALCULATE SOLUTION FOR IMPORT NON EXEC, the owning rank will do that
    for (std::size_t i = 0; i < n_nodes_set; i++) 
        if (((int*)node_type_dat->data)[i] == NORMAL || ((int*)node_type_dat->data)[i] == OPEN) 
            neq++;

    ID = new int[n_nodes_inc_halo];
   
    LM = new int*[n_elements_inc_halo]; /*allocate location matrix, n_elements*4 */
    for (int e = 0; e < n_elements_inc_halo; e++) 
        LM[e] = new int[4];

    NX = new double**[n_elements_inc_halo]; /*allocate NX matrix*/
    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        NX[e] = new double*[4];
        for (int a=0;a<4;a++) 
            NX[e][a] = new double[3];
    }
    
    d = new double[neq];        /*solution array*/
    for (int n=0;n<neq;n++) 
        d[n]=0;                 /*initial guess*/

    detJ = new double[n_elements_inc_halo];

    /*set quadrature points*/
    l[0]=-sqrt(1.0/3.0); 
    l[1]=sqrt(1.0/3.0);
    W[0]=1; 
    W[1]=1;
    n_int = 2;

    initPetscStructures();

    initID(node_type_dat);

    for (int e = 0; e < n_elements_inc_halo; e++) /*now set up the LM matrix*/
    {
        for (int a=0;a<4;a++)    /*tetrahedra*/ 
        {
            int node_idx = cell_to_nodes_map->map[e * cell_to_nodes_map->dim + a];
            LM[e][a] = ID[node_idx]; 
        }
    }

    computeNX(node_pos_dat, cell_to_nodes_map);

    summarize(std::cout);

    preAssembly(cell_to_nodes_map, node_bnd_pot_dat);

    if (print_petsc) 
    {
        opp_printf("FESolver::FESolver", "This is KMat ****************************************");
        MatView(Kmat, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver::FESolver", "Above is KMat ****************************************");
    }
    if (print_petsc) 
    {
        opp_printf("FESolver::FESolver", "This is F0vec ****************************************");   
        VecView(F0vec, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver::FESolver", "Above is F0vec ****************************************");
    }

    sanityCheck();
}

//*************************************************************************************************
FESolver::~FESolver() 
{
    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        for (int a=0;a<4;a++) 
            delete[] NX[e][a];
        delete[] NX[e];
        delete[] LM[e];
    }

    delete[] LM;
    delete[] NX;
    delete[] ID;
    delete[] d;
    delete[] detJ;

    /* Petsc */
    delete[] vecCol;
    delete[] matCol;
    for (int i=0; i<neq; i++) { delete[] matIndex[i]; delete[] matIndexValues[i]; }
    delete[] matIndex;
    delete[] matIndexValues;
    
    KSPDestroy(&ksp);
    VecDestroy(&Xvec);
    VecDestroy(&Bvec);
    VecDestroy(&F0vec);
    VecDestroy(&F1vec);
    VecDestroy(&Dvec);
    VecDestroy(&Yvec);
    MatDestroy(&Jmat);
    MatDestroy(&Kmat);
}

//*************************************************************************************************
void FESolver::initID(oppic_dat node_type_dat) // relation from node indices to equation indices
{
    int P=0;
    for (int n = 0; n < n_nodes_inc_halo; n++)
    {
        if (n < n_nodes_set && 
            (((int*)node_type_dat->data)[n] == NORMAL || ((int*)node_type_dat->data)[n] == OPEN))
        {
            ID[n]=P;
            P++;
        }
        else 
        {
            ID[n]=INT_MIN;    /* dirichlet node or for halo region, until updated by MPI comm */
        }
    }

#ifdef ENABLE_MPI
    // send all the ID values of the export non-exec halos
    // (send the values by adding by own_start and making the values negative)
    // receive all ID values of the import non exec halos and 
    // assign the values to the correct indices of ID array

    halo_list exp_nonexec_list = OP_export_nonexec_list[node_type_dat->set->index];
    halo_list imp_nonexec_list = OP_import_nonexec_list[node_type_dat->set->index];

    int* send_buffer = new int[exp_nonexec_list->size];
    std::vector<MPI_Request> send_req(exp_nonexec_list->ranks_size);

    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) 
    {
        for (int j = 0; j < exp_nonexec_list->sizes[i]; j++) 
        {
            int index = exp_nonexec_list->disps[i] + j;           
            send_buffer[index] = ID[exp_nonexec_list->list[index]];

            if (send_buffer[index] != INT_MIN)
                send_buffer[index] = -1 * (send_buffer[index] + own_start);

            // opp_printf("SEND", "idx=%d val=%d to rank %d", exp_nonexec_list->list[index], 
            //     send_buffer[index], exp_nonexec_list->ranks[i]);
        }

        MPI_Isend(&send_buffer[exp_nonexec_list->disps[i]], exp_nonexec_list->sizes[i], MPI_INT,
                exp_nonexec_list->ranks[i], 11, OP_MPI_WORLD, &(send_req[i]));
    }

    int nonexec_init = node_type_dat->set->size + node_type_dat->set->exec_size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    {
        // opp_printf("SEND", "count=%d from rank %d", imp_nonexec_list->sizes[i], imp_nonexec_list->ranks[i]);

        MPI_Recv(&(ID[nonexec_init + imp_nonexec_list->disps[i]]), imp_nonexec_list->sizes[i], MPI_INT,
            imp_nonexec_list->ranks[i], 11, OP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(exp_nonexec_list->ranks_size, send_req.data(), MPI_STATUSES_IGNORE);

    send_req.clear();
    delete[] send_buffer;
#endif
}

//*************************************************************************************************
void FESolver::initPetscStructures()
{
    MatCreate(PETSC_COMM_WORLD, &Jmat);
    MatCreate(PETSC_COMM_WORLD, &Kmat);

    MatSetSizes(Jmat, neq, neq, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetSizes(Kmat, neq, neq, PETSC_DETERMINE, PETSC_DETERMINE);

#ifndef ENABLE_MPI
    MatSetType(Jmat, MATSEQAIJ);
    MatSetType(Kmat, MATSEQAIJ);

    VecCreateSeq(PETSC_COMM_WORLD, neq, &Bvec);
#else
    MatSetType(Jmat, MATMPIAIJ);
    MatSetType(Kmat, MATMPIAIJ);

    VecCreateMPI(PETSC_COMM_WORLD, neq, PETSC_DETERMINE, &Bvec);
#endif

    MatSetFromOptions(Jmat);
    MatSetFromOptions(Kmat);
    
    MatSetUp(Jmat);
    MatSetUp(Kmat);

    VecSetFromOptions(Bvec);
    
    VecDuplicate(Bvec, &Xvec);
    VecDuplicate(Bvec, &F0vec);
    VecDuplicate(Bvec, &F1vec);
    VecDuplicate(Bvec, &Dvec);
    VecDuplicate(Bvec, &Yvec);

    VecSet(Bvec, 0.0);
    VecSet(F0vec, 0.0);
    VecSet(F1vec, 0.0);
    VecSet(Dvec, 0.0);
    VecSet(Yvec, 0.0);
    VecAssemblyBegin(Bvec); VecAssemblyEnd(Bvec);
    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);
    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);
    VecAssemblyBegin(Dvec); VecAssemblyEnd(Dvec);
    VecAssemblyBegin(Yvec); VecAssemblyEnd(Yvec);

    VecGetOwnershipRange(Bvec, &own_start, &own_end);
    VecGetSize(Bvec, &global_neq);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, Jmat, Jmat);
    // KSPSetTolerances(ksp, 1.e-2 / (neq * neq), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); 
    // KSPSetTolerances(ksp, 1.e-100 / (neq * neq), 1.e-100, PETSC_DEFAULT, PETSC_DEFAULT); 
    KSPSetTolerances(ksp, 1.e-100, 1.e-100, PETSC_DEFAULT , PETSC_DEFAULT); 
    KSPSetFromOptions(ksp); 

    //PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

    vecCol         = new int[neq];
    matCol         = new int[neq];
    matIndex       = new int*[neq];
    matIndexValues = new double*[neq];

    for (int i = 0; i < neq; i++)
    {
        vecCol[i] = i + own_start;
    }
}

//*************************************************************************************************
/* preassembles the K matrix and "F0 vector */
void FESolver::preAssembly(oppic_map cell_to_nodes_map, oppic_dat node_bnd_pot) 
{

    opp_printf("FESolver", "preAssembly");

    double **K = new double*[neq]; // K will be (neq x global_neq) matrix
    for (int i=0;i<neq;i++)
    {
        K[i] = new double[global_neq];
        
        for (int j=0;j<global_neq;j++) 
            K[i][j] = 0;
    }

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        double ke[4][4];

        for (int a=0;a<4;a++)
            for (int b=0;b<4;b++) 
            {
                ke[a][b] = 0;    /*reset*/

                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            double nax[3],nbx[3];

                            getNax(nax,e,a);
                            getNax(nbx,e,b);

                            /*dot product*/
                            double dot=0;
                            for (int d=0;d<3;d++) dot+=nax[d]*nbx[d];
                            ke[a][b] += dot*detJ[e]*W[i]*W[j]*W[k];
                        }
            }

        addKe(K, e,ke);
   
        double fe[4]; /*force vector*/

        for (int a=0;a<4;a++) 
        {  
            double fh=0;    /*second term int(na*h), always zero since support only h=0*/
            double fg = 0;  /*third term, -sum(kab*qb)*/

            for (int b=0;b<4;b++) 
            {
                // This can reach import halos
                int node_idx = cell_to_nodes_map->map[e * cell_to_nodes_map->dim + b]; 
                double gb = ((double*)node_bnd_pot->data)[node_idx];
                fg-=ke[a][b]*gb;
            }

            fe[a] = fh + fg; /*combine*/
        }

        addFe(&F0vec, e,fe);

    }  /*end of element*/

    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);

    initialzeMatrix(K);

    double scalar1 = 1 / SCALLING;
    double scalar2 = 1 / (SCALLING * SCALLING); 

    VecScale(F0vec, scalar2);   // downscalling since NX is scalled to avoid precision issues
    MatScale(Kmat, scalar2);    // downscalling since NX is scalled to avoid precision issues

    for (int i=0;i<neq;i++) 
        delete[] K[i];
    delete[] K;
}

//*************************************************************************************************
void FESolver::buildF1Vector(double *ion_den) 
{
    double *f = new double[neq];

    for (int n = 0; n < n_nodes_inc_halo; n++) 
    {
        int A = ID[n];
        if (A<0) continue;    /*skip known nodes*/
        f[A] = (QE/EPS0)*((ion_den[n]) + n0*exp((d[A]-phi0)/kTe));
    }

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        double fe[4];
        for (int a=0;a<4;a++) 
        {
            /*first term is int(na*f), set to zero for now*/
            double ff=0;
            int A = LM[e][a];

            // IDEA : F1vec will be enriched only for the own range (owner compute),
            // hence import halo region calculation is not required here. is it correct?
            if (A>=0)    /*if unknown node*/ 
            {
                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += f[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }

                ff*=(1.0/8.0);    /*change of limits*/
                fe[a] = ff;
            }
        }

        addFe(&F1vec, e,fe);
    }

    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);

    delete[] f;
}

//*************************************************************************************************
void FESolver::buildJmatrix() 
{ 

    double *fp_term = new double[neq];
    double *FP = new double[neq];

    for (int n=0;n<neq;n++) 
        FP[n] = 0;

    for (int n=0;n<neq;n++) 
        fp_term[n] = -QE/EPS0*n0*exp((d[n]-phi0)/kTe)*(1/kTe);

    MatCopy(Kmat, Jmat, DIFFERENT_NONZERO_PATTERN);

    double fe[4];  /*build fprime vector*/

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        for (int a=0;a<4;a++) 
        {
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/ 
            {
                for (int k=0;k<n_int;k++) /*perform quadrature*/
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += fp_term[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }

                ff*=(1.0/8.0);    /*change of limits*/
            }

            fe[a] = ff;
        }

        /*assembly*/
        for (int a=0;a<4;a++)    /*tetrahedra*/ 
        {
            int P = LM[e][a];
            if (P<0) continue;    /*skip g nodes*/

            FP[P] += fe[a];
        }
    }

    for (int u=0;u<neq;u++) /*subtract diagonal term*/
    {
        MatSetValue(Jmat, (u + own_start), (u + own_start), (-FP[u]), ADD_VALUES); // J[u][u]-=FP[u];
    }

    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    delete[] fp_term;
    delete[] FP;
}


//*************************************************************************************************
/*wrapper for solving the non-linear Poisson's equation*/
void FESolver::computePhi(oppic_arg arg0, oppic_arg arg1, oppic_arg arg2) 
{ TRACE_ME;

    if (OP_DEBUG) opp_printf("FESolver", "computePhi START");

    opp_profiler->start("ComputePhi");

    int nargs = 3;
    oppic_arg args[nargs];
    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    args[1].idx = 2; // HACK to forcefully make halos to download

    int set_size = opp_mpi_halo_exchanges(arg1.dat->set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    if (false) // incase if we want to print current ranks node_charge_density including halos
    {
        std::string f = std::string("FESolverComputePhi_") + std::to_string(OPP_rank);
        opp_print_dat_to_txtfile(args[1].dat  , f.c_str(), "node_charge_density.dat");
    }

    double *ion_den = (double*)arg1.dat->data;
    double *node_potential = (double*)arg0.dat->data;
    double *node_bnd_pot = (double*)arg2.dat->data;

    if (fesolver_method == Method::NonLinear) 
    {
        nonLinearSolve(ion_den);
    } 
    else if (fesolver_method == Method::Petsc) 
    {
        if (linearSolve(ion_den) == false)
        {
            char* str_reason = "\0";
            KSPGetConvergedReasonString(ksp, (const char**)(&str_reason));
            std::cerr << "linearSolve Petsc failed to converge : " << str_reason <<
                "; run with -ksp_converged_reason" << std::endl;           
        }
    }
    else
    {
        opp_printf("FESolver", "Error... Solver type not implemented");
    }

    /* combine d (solution from linear solver) and boundary potential to get node_potential */
    for (int n = 0;n < n_nodes_set; n++) // owner compute, hence can have only upto node set size
    {
        node_potential[n] = node_bnd_pot[n]; /*zero on non-g nodes*/

        int A=ID[n];
        if (A>=0)           /*is this a non-boundary node?*/
            node_potential[n] += d[A];
    }

    opp_mpi_set_dirtybit(nargs, args);

    opp_profiler->end("ComputePhi");

    if (OP_DEBUG) opp_printf("FESolver", "computePhi DONE");
}

//*************************************************************************************************
void FESolver::nonLinearSolve(double *ion_den) 
{ TRACE_ME;

    int it = 0;

    if (OP_DEBUG) opp_printf("FESolver", "nonLinearSolve START");
    bool converged = false;

    for (; it<10; it++) 
    {
        if (converged = linearSolve(ion_den))
            break;
    }

    if (!converged) 
    {
        char* str_reason = "\0";
        KSPGetConvergedReasonString(ksp, (const char**)(&str_reason));
        std::cerr << "nonLinearSolve Petsc failed to converge : " << str_reason << 
            "; run with -ksp_converged_reason"<< std::endl;
    }

    if (OP_DEBUG) opp_printf("FESolver", "nonLinearSolve DONE it=%d", it);
}

//*************************************************************************************************
bool FESolver::linearSolve(double *ion_den) 
{ TRACE_ME;

    if (OP_DEBUG) opp_printf("FESolver", "linearSolve START");

    buildF1Vector(ion_den);

    if (print_petsc) 
    {
        opp_printf("FESolver", "This is F1vec ***************************************************");
        VecView(F1vec, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver", "Above is F1vec ***************************************************");
    }

    MatMult(Kmat, Dvec, Bvec);  // B = K * D
    
    VecAXPY(Bvec, -1.0, F0vec); // B = B - F0

// TODO_IMM : Can build F1 vec here by waiting for halo exch complete. This will overlap some work

    VecAXPY(Bvec, -1.0, F1vec); // B = B - F1

    buildJmatrix();             // Build J using K and D

    if (print_petsc) 
    {
        opp_printf("FESolver", "This is J ***************************************************");
        MatView(Jmat, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver", "Above is J ***************************************************");
    }

    if (print_petsc) 
    {
        opp_printf("FESolver", "This is B ***************************************************");
        VecView(Bvec, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver", "Above is B ***************************************************");
    }

    KSPSolve(ksp, Bvec, Yvec); // Jmat * Yvec = Bvec (Yvec = solution)

    if (print_petsc) 
    {
        opp_printf("FESolver", "This is Y ***************************************************");
        VecView(Yvec, PETSC_VIEWER_STDOUT_WORLD);
        opp_printf("FESolver", "Above is Y ***************************************************");
    }

    VecAXPY(Dvec, -1.0, Yvec); // D = D - Y

    VecGetValues(Dvec, neq, vecCol, d); // For the calculation at computePhi()

    KSPGetConvergedReason(ksp, &reason);

    if (OP_DEBUG) opp_printf("linearSolve", "solve DONE");

    return (reason >= 0);
}

//*************************************************************************************************
void FESolver::initialzeMatrix(double **p_A)
{
    if (OP_DEBUG) opp_printf("FESolver", "initialzeMatrix START");

    if (!matIndexCreated)
    {
        int diag_max_fields = 0, off_diag_max_fields = 0;
        
        for (int i=own_start; i<own_end; i++) // iterate own rows
        {
            std::vector<int> tempVec;
            int local_diag_max_fields = 0, local_off_diag_max_fields = 0;
            
            for (int j=0; j<global_neq; j++) // iterate all columns
            {    
                // significant ones and all diagonals         
                if ((std::abs(p_A[i-own_start][j]) > 1e-12) || (i == j)) 
                {
                    tempVec.push_back(j); 

                    if (j >= own_start && j < own_end)
                        local_diag_max_fields++;
                    else
                        local_off_diag_max_fields++;
                }         
            }

            matCol[i-own_start] = tempVec.size();
            matIndex[i-own_start] = new int[tempVec.size()];
            matIndexValues[i-own_start] = new double[tempVec.size()];

            std::copy(tempVec.begin(), tempVec.end(), matIndex[i-own_start]);

            diag_max_fields = (diag_max_fields > local_diag_max_fields) ? 
                                    diag_max_fields : local_diag_max_fields; 

            off_diag_max_fields = (off_diag_max_fields > local_off_diag_max_fields) ? 
                                    off_diag_max_fields : local_off_diag_max_fields;  
        }

#ifndef ENABLE_MPI
        MatSeqAIJSetPreallocation(Kmat, diag_max_fields, nullptr);
        MatSeqAIJSetPreallocation(Jmat, diag_max_fields, nullptr);
#else
        MatMPIAIJSetPreallocation(Kmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
        MatMPIAIJSetPreallocation(Jmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
#endif

        if (OP_DEBUG) opp_printf("FESolver", "initialzeMatrix diag_max_fields=%d off_diag_max_fields=%d", 
            diag_max_fields, off_diag_max_fields);
        matIndexCreated = true;
    }

    for (int i=own_start; i<own_end; i++)
    {
        for (int j=0; j<matCol[i-own_start]; j++) 
            matIndexValues[i-own_start][j] = p_A[i-own_start][matIndex[i-own_start][j]];

        MatSetValues(Kmat, 1, &i, matCol[i-own_start], matIndex[i-own_start], 
            matIndexValues[i-own_start], INSERT_VALUES); 
    }

    MatAssemblyBegin(Kmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Kmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    if (OP_DEBUG) opp_printf("FESolver", "initialzeMatrix END");
}

//*************************************************************************************************

// UTIL FUNCTIONS

//*************************************************************************************************
void FESolver::summarize(std::ostream &out) 
{
    opp_printf("FESolver", "FE SOLVER INFORMATION");
    opp_printf("FESolver", "---------------------");
    opp_printf("FESolver", "own neq is %d, the global neq is %d, start_eq is %d, end_eq %d -------", 
        neq, global_neq, own_start, own_end);
    opp_printf("FESolver", "---------------------");
}

void FESolver::enrich_cell_shape_deriv(oppic_dat cell_shape_deriv)
{
    if (OP_DEBUG) opp_printf("FESolver", "enrich_cell_shape_deriv START");

    opp_profiler->start("EnrichCellShapeDeriv");

    // copying only up to set size, hence import exec halos will be zero
    for (int cellID = 0; cellID < n_elements_inc_halo; cellID++)
    {
        for (int nodeCon = 0; nodeCon < N_PER_C; nodeCon++)
        {
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 0 ] = 
                (NX[cellID][nodeCon][0] / SCALLING);
            
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 1 ] = 
                (NX[cellID][nodeCon][1] / SCALLING);
            
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 2 ] = 
                (NX[cellID][nodeCon][2] / SCALLING);
        }
    }

#ifdef ENABLE_MPI
    cell_shape_deriv->dirtybit = 1;
    oppic_arg arg0 = opp_get_arg(cell_shape_deriv, OP_RW);
    arg0.idx = 2; // HACK
    opp_mpi_halo_exchanges(cell_shape_deriv->set, 1, &arg0);  
    opp_mpi_halo_wait_all(1, &arg0);
#endif

    opp_profiler->end("EnrichCellShapeDeriv");

    if (OP_DEBUG) opp_printf("FESolver", "enrich_cell_shape_deriv END");
}

/*adds contributions from element stiffness matrix*/
void FESolver::addKe(double** K, int e, double ke[4][4]) // BUG : K is not created correctly
{ 

    if (neq <= 0) return;

    for (int a=0;a<4;a++)         /*tetrahedra*/
    {
        for (int b=0;b<4;b++) 
        {
            int P = LM[e][a];
            int Q = LM[e][b];
            if (P<0) continue;    /* skip g nodes or not in current ranks own row range */

            if (Q<0)
            {
                if (Q == INT_MIN) continue; // skip g nodes
                else Q = -1 * Q;            // on a different MPI rank, get its correct column
            }
            else 
                Q += own_start;             // on the same MPI rank, get its correct column

            K[P][Q] += ke[a][b];
        }
    }
}

/*adds contributions from element force vector to a global F vector*/
void FESolver::addFe(Vec *Fvec, int e, double fe[4]) 
{ 

    if (neq <= 0) return;
    
    for (int a=0;a<4;a++)    /*tetrahedra*/ 
    {
        int P = LM[e][a];
        if (P<0) continue;    /* skip g nodes or on a different rank (owner compute)*/
        
        VecSetValue(*Fvec, P + own_start, fe[a], ADD_VALUES); // F[P] += fe[a];
    }
}

/*evaluates shape function a at position (xi,eta,zeta)*/
double FESolver::evalNa(int a, double xi, double eta, double zeta) 
{
    switch(a) 
    {
        case 0: return xi;
        case 1: return eta;
        case 2: return zeta;
        case 3: return 1-xi-eta-zeta;       
    }

    return 0;    /*shouldn't happen*/
}

/*returns derivative of N[a] at some logical point since we are using linear elements, 
these are constant in each element*/
void FESolver::getNax(double nx[3], int e, int a) 
{
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

/*computes derivatives of the shape functions for all elements constants since 
using linear elements*/
void FESolver::computeNX(oppic_dat node_pos, oppic_map cell_to_nodes_map) 
{ 

    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {       
        /*node indices*/
        int* map0idx = &((int*)cell_to_nodes_map->map)[e * cell_to_nodes_map->dim]; 

        double x[4][3];
        for (int a=0;a<4;a++) 
        {
            /*node positions*/
            double *pos = &((double*)node_pos->data)[map0idx[a] * node_pos->dim]; 
            for (int d=0;d<3;d++) 
                x[a][d] = pos[d];
        }

        double x_xi[3][3];

        for (int i=0;i<3;i++)    /*x/y/z*/
            for (int j=0;j<3;j++) /*xi/eta/zeta*/ 
            {
                x_xi[i][j] = 0;
                for (int a=0; a<4; a++)    /*tet node*/
                    x_xi[i][j] += na_xi[a][j]*x[a][i];
            }

        detJ[e] = det3(x_xi);

        double xi_x[3][3];
        inverse(x_xi,xi_x);
      
        for (int a=0;a<4;a++) /*evaluate na_x*/
            for (int d=0;d<3;d++) 
            {
                NX[e][a][d]=0;
                for (int k=0;k<3;k++)
                    // scalled to avoid precision issues
                    NX[e][a][d] += (na_xi[a][k]*xi_x[k][d] * SCALLING); 
            }
    }
}

//*************************************************************************************************
// Sanity check after assembling matrices
void FESolver::sanityCheck()
{
    int j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end;

    MatGetSize(Jmat, &j_gm, &j_gn);
    MatGetSize(Kmat, &k_gm, &k_gn);

    MatGetOwnershipRange(Jmat, &j_row_start, &j_row_end);
    MatGetOwnershipRange(Kmat, &k_row_start, &k_row_end);

    if (is_neq(j_gm, j_gn) || is_neq(k_gm, k_gn) || is_neq(j_gm, global_neq) ||
        is_neq(j_row_start, k_row_start) || is_neq(j_row_end, k_row_end) || 
        is_neq(j_row_start, own_start) || is_neq(j_row_end, own_end))
    {
        opp_printf("FESolver::FESolver", "Error... Matrix vector sizes issue");
    }      
    
    if (OP_DEBUG) 
        opp_printf("FESolver::FESolver", "j_gm %d, j_gn %d, k_gm %d, k_gn %d, j_row_start %d, \
            j_row_end %d, k_row_start %d, k_row_end %d", 
            j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end);
}

//*************************************************************************************************