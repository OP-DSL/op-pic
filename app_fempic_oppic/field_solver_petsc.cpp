// TODO : To be converted to API with kernels

#include "fempic.h"

//*************************************************************************************************
/*FESolver*/
FESolver::FESolver(Volume &volume, int argc, char **argv)
    : volume(volume)
{
    cout<<"FESolver PETSC RUNNING : ************************"<<endl;
    
    int dummy = 0;
    //PetscInitialize(&dummy, NULL, (char *)0, "petsc FESolver");
    PetscInitialize(&argc, &argv, PETSC_NULL, "petsc FESolver PP");

    /*count number of unknowns*/
    neq     = 0;
    phi0    = 0;
    n0      = PLASMA_DEN;
    kTe     = 2;

    /*OPEN nodes are "h" nodes*/
    for (size_t i=0; i<volume.nodes.size(); i++)
        if (volume.nodes[i].type == NORMAL || volume.nodes[i].type == OPEN) //except INLET and SPHERE
            neq++;
    cout<<"FESolver : There are "<<neq<<" unknowns"<<endl;
    
    /*allocate neq*neq K matrix*/
    p_K = new double*[neq];
    for (int i=0;i<neq;i++) 
        p_K[i] = new double[neq];
    cout<<"FESolver : Allocated "<<neq<<"x"<<neq<<" stiffness matrix"<<endl;

    /*allocate neq*neq J matrix*/
    p_J = new double*[neq];
    for (int i=0; i<neq; i++) 
        p_J[i] = new double[neq];
    cout<<"FESolver : Allocated "<<neq<<"x"<<neq<<" Jacobian matrix"<<endl;

    /*allocate F0 and F1 vectors*/
    F0 = new double[neq];
    F1 = new double[neq];
    cout<<"FESolver : Allocated two "<<neq<<"x1 force vectors"<<endl;

    n_nodes = volume.nodes.size();
    n_cells = volume.cells.size();

    /*allocate ID vector*/
    ID = new int[n_nodes];
    cout<<"FESolver : Allocated "<<n_nodes<<"x1 ID vector"<<endl;

    /*allocate location matrix, n_cells*4 */
    LM = new int*[n_cells];
    for (int e=0; e<n_cells; e++) 
        LM[e] = new int[4];
    cout<<"FESolver : Allocated "<<n_cells<<"x4 location matrix"<<endl;

    /*allocate NX matrix*/
    NX = new double**[n_cells];
    for (int e=0; e<n_cells; e++) 
    {    
        NX[e] = new double*[4];
        for (int a=0; a<4; a++) 
            NX[e][a] = new double[3];
    }
    cout << "FESolver : Allocated " << n_cells << "x4x3 NX matrix" << endl;

    /*solution array*/
    p_d = new double[neq];
    for (int n=0;n<neq;n++) 
        p_d[n]=0;    /*initial guess*/
    
    detJ = new double[n_cells];
    cout<<"FESolver : Allocated "<<n_cells<<"x1 detJ vector"<<endl;

    /*set up the ID array
    note valid values are 0 to neq-1 and -1 indicates "g" node*/
    int P=0;
    for (int n=0; n<n_nodes; n++)
    {
        if (volume.nodes[n].type==NORMAL ||  volume.nodes[n].type==OPEN) 
        {
            ID[n]=P;
            P++;
        }
        else
        {
            ID[n]=-1;    /*dirichlet node*/
        }
    }

    /*now set up the LM matrix*/
    for (int e=0; e<n_cells; e++)
    {
        for (int a=0; a<4; a++)    /*tetrahedra*/ // Iterate over the attached 4 nodes of the cell, (volume.cells[e].con[a])
        {
            LM[e][a] = ID[volume.cells[e].con[a]];
        }
    }
    cout<<"FESolver : Built ID and LM matrix"<<endl;
    
    /*set quadrature points*/
    p_l[0] = -sqrt(1.0/3.0); 
    p_l[1] = sqrt(1.0/3.0);

    p_W[0] = 1; 
    p_W[1] = 1;

    n_int = 2;

    cout<<"FESolver : computeNX"<<endl;
    computeNX();

    /* Petsc */
    vecCol           = new int[neq];
    matCol           = new int[neq];
    matIndex         = new int*[neq];
    matIndexValues   = new double*[neq];

    for (int i = 0; i < neq; i++)
    {
        vecCol[i] = i;
        matCol[i] = 0;
    }

    cout << "FESolver : Petsc Creating Matrix " << endl;

    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, neq, neq);
    MatSetFromOptions(A);
    MatSetUp(A);

    cout << "FESolver : Petsc Creating Vector " << endl;

    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetSizes(b, PETSC_DECIDE, neq);
    VecSetFromOptions(b);
    VecDuplicate(b, &x);

    cout << "FESolver : Petsc Creating KSP " << endl;

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    // KSPSetTolerances(ksp, 1.e-2 / (neq * neq), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); 
    KSPSetTolerances(ksp, 1.e-100 / (neq * neq), 1.e-100, PETSC_DEFAULT, 100); 
    KSPSetFromOptions(ksp);

    cout << "FESolver : Petsc Done " << endl;
}

//*************************************************************************************************
/*~FESolver, frees memory*/
FESolver::~FESolver()
{
    for (int i=0;i<neq;i++) {delete[] p_K[i]; delete[] p_J[i];}
    for (int e=0;e<n_cells;e++) delete[] LM[e];
    for (int e=0;e<n_cells;e++) 
    {    
        for (int a=0;a<4;a++) delete[] NX[e][a];
        delete NX[e];
    }

    delete[] p_K;
    delete[] p_J;
    delete[] LM;
    delete[] F0;
    delete[] F1;
    delete[] NX;
    delete[] ID;
    delete[] p_d;
    delete[] detJ;

    /* Petsc */
    delete[] vecCol;
    delete[] matCol;
    for (int i=0; i<neq; i++) { delete[] matIndex[i]; delete[] matIndexValues[i];}
    delete[] matIndex;
    delete[] matIndexValues;
    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

    PetscFinalize();
}

//*************************************************************************************************
/*clears K and F*/
void FESolver::startAssembly()
{
    for (int i=0;i<neq;i++)
        for (int j=0;j<neq;j++) 
            p_K[i][j] = 0;

    for (int i=0;i<neq;i++) {
        F0[i]=0;
        F1[i]=0;
    }
}

//*************************************************************************************************
/*adds contributions from cell stiffness matrix*/
void FESolver::addKe(int e, double ke[4][4])
{
    for (int a=0; a<4; a++)    /*tetrahedra*/
    {
        for (int b=0; b<4; b++)
        {
            int P = LM[e][a];
            int Q = LM[e][b];
            if (P<0 || Q<0) continue;    /*skip g nodes*/

            p_K[P][Q] += ke[a][b];
        }
    }
}

//*************************************************************************************************
/*adds contributions from cell force vector to a global F vector*/
void FESolver::addFe(double *F, int e, double fe[4])
{
    for (int a=0; a<4; a++)    /*tetrahedra*/
    {
        int P = LM[e][a];
        if (P<0) continue;    /*skip g nodes*/

        F[P] += fe[a];
    }
}

//*************************************************************************************************
/*evaluates shape function a at position (xi,eta,zeta)*/
double FESolver::evalNa(int a, double xi, double eta, double zeta)
{
    switch(a)
    {
    case 0: return xi; break;
    case 1: return eta; break;
    case 2: return zeta; break;
    case 3: return 1-xi-eta-zeta; break;
    default: return 0;    /*shouldn't happen*/
    }
}

//*************************************************************************************************
/*returns derivative of N[a] at some logical point
since we are using linear cells, these are constant in each cell*/
void FESolver::getNax(double nx[3], int e, int a, double xi, double eta, double zeta)
{
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

//*************************************************************************************************
/*computes derivatives of the shape functions for all cells
constants since using linear cells*/
void FESolver::computeNX()
{
    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};
    
    for (int e=0; e<n_cells; e++)
    {
        /*node positions*/
        Tetra &tet = volume.cells[e];

        double x[4][3];
        for (int a=0;a<4;a++)
        {
            double *pos = volume.nodes[tet.con[a]].pos;
            for (int d=0;d<3;d++) 
                x[a][d] = pos[d];    /*copy*/
        }

        /*compute x_xi matrix*/
        double x_xi[3][3];

        for (int i=0;i<3;i++)    /*x/y/z*/
            for (int j=0;j<3;j++) /*xi/eta/zeta*/
            {
                x_xi[i][j] = 0;
                for (int a=0; a<4; a++)    /*tet node*/
                    x_xi[i][j] += na_xi[a][j]*x[a][i];
            }

        /*save det(x_xi)*/
        detJ[e] = det3(x_xi);

        /*compute matrix inverse*/
        double xi_x[3][3];
        inverse(x_xi,xi_x);

        /*evaluate na_x*/
        for (int a=0;a<4;a++)
        {
            for (int d=0;d<3;d++)
            {
                NX[e][a][d]=0;
                for (int k=0;k<3;k++)
                    NX[e][a][d]+=na_xi[a][k]*xi_x[k][d];    
            }
        }
    }
}

//*************************************************************************************************
/*compute inverse of a 3x3 matrix using the adjugate method*/
void FESolver::inverse(double M[3][3], double V[3][3])
{
    double a=M[0][0];
    double b=M[0][1];
    double c=M[0][2];
    double d=M[1][0];
    double e=M[1][1];
    double f=M[1][2];
    double g=M[2][0];
    double h=M[2][1];
    double i=M[2][2];

    V[0][0]=(e*i-f*h);
    V[1][0]=-(d*i-f*g);
    V[2][0]=(d*h-e*g);
    V[0][1]=-(b*i-c*h);
    V[1][1]=(a*i-c*g);
    V[2][1]=-(a*h-b*g);
    V[0][2]=(b*f-c*e);
    V[1][2]=-(a*f-c*d);
    V[2][2]=(a*e-b*d);
    double det = a*V[0][0]+b*V[1][0]+c*V[2][0];

    double idet=0;
    if (fabs(det)<1e-12) {
        cerr<<"FESolver : Matrix is not invertible, setting to [0]!"<<endl;}
    else idet=1/det;
    
    /*1/det*/
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            V[i][j]*=idet;
}

//*************************************************************************************************
/*Newton Rhapson solver, input is the ion density, writing to d*/
void FESolver::solveNonLinear(const double *ion_den)
{
    const double tol = 1e-2;

    /*allocate memory for y*/
    double *y = new double[neq];
    double *G = new double[neq];
    
    /*clear y values*/
    for (int i=0; i<neq; i++) 
        y[i]=0;
    
    bool converged = false;
    double L2;
    for (int it=0; it<10; it++)
    {
        /*builds the "ff" part of the force vector*/
        buildF1Vector(ion_den);         // F0 vector is already built and building F1 vector here

        /*form G=K*d-F*/
        /* CAN */ matVecMultiply(G,p_K,p_d, neq);    //G=K*d, writing to G // K is the global stiffness matrix
        /* CAN */ vecVecSubtract(G,G,F0, neq);       //G=G-F0 giving us G=G-F0
        /* CAN */ vecVecSubtract(G,G,F1, neq);       //G=G-F1 giving us G=G-F1

        buildJmatrix();

        /* CAN PCSOR or even with KSP*/ 
        // solveLinear(p_J,y,G); /*simple Gauss-Seidel solver for J*y=G, writing to y*/
        solveLinear_petsc(p_J,y,G);
        
        /*now that we have y, update solution */
        for (int n=0; n<neq; n++) 
            p_d[n] -= y[n];

        /*compute residue*/
        double sum=0;
        for (int u=0; u<neq; u++)
        {
            sum += y[u]*y[u];
        }
        L2 = sqrt(sum)/neq;

        if (L2 < tol) 
        {
            cout<<" NR converged in "<<it+1<<" iterations with L2="<<setprecision(3)<<L2<<endl;
            converged=true;
            break;
        }
    }

    delete[] y;
    delete[] G;

    if (!converged) 
    {
        cerr<<"FESolver : NR failed to converge, L2 = "<<L2<<endl;
        exit(-10);
    }            
}

//*************************************************************************************************
/*builds J matrix for NR solver*/
void FESolver::buildJmatrix()
{
    /*first compute exponential term*/
    double *fp_term = new double[neq];
    double *FP = new double[neq];

    for (int n=0;n<neq;n++) 
        FP[n] = 0;

    for (int n=0;n<neq;n++)
    {
        fp_term[n] = -QE/EPS0*n0*exp((p_d[n]-phi0)/kTe)*(1/kTe);
    }

    /*now set J=K*/
    for (int i=0; i<neq; i++)
    {
        for (int j=0; j<neq; j++)
        {
            p_J[i][j] = p_K[i][j];
        }
    }

    /*build fprime vector*/
    double fe[4];

    for (int e=0; e<n_cells; e++)
    {
        for (int a=0; a<4; a++)
        {    
            double ff=0;                
            int A = LM[e][a];    // Gives the equation number with respect to cell e and its a th node
            if (A>=0)    /*if unknown node*/
            {
                /*perform quadrature*/
                for (int k=0; k<n_int; k++)
                    for (int j=0; j<n_int; j++)
                        for (int i=0; i<n_int; i++)
                        {
                            /*change of limits*/
                            double xi = 0.5*(p_l[i]+1);
                            double eta = 0.5*(p_l[j]+1);
                            double zeta = 0.5*(p_l[k]+1);

                            double Na = evalNa(a,xi,eta,zeta);

                            ff += fp_term[A]*Na*detJ[e]*p_W[i]*p_W[j]*p_W[k];
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


    /*subtract diagonal term*/
    for (int u=0;u<neq;u++)
    {
        p_J[u][u]-=FP[u];
    }

    delete[] fp_term;
    delete[] FP;
}

//*************************************************************************************************
/*preassembles the K matrix and "h" and "g" parts of the force vector*/
void FESolver::preAssembly(double *boundary_potential)
{
    /*loop over cells*/
    for (int e=0;e<n_cells;e++)
    {
        Tetra &tet = volume.cells[e];        
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

                            double xi = 0.5*(p_l[i]+1);    /*not used!*/
                            double eta = 0.5*(p_l[j]+1);
                            double zeta = 0.5*(p_l[k]+1);
                            getNax(nax,e,a,xi,eta,zeta);
                            getNax(nbx,e,b,xi,eta,zeta);
                            
                            /*dot product*/
                            double dot=0;
                            for (int d=0;d<3;d++) dot+=nax[d]*nbx[d];
                            ke[a][b] += dot*detJ[e]*p_W[i]*p_W[j]*p_W[k];
                        }
            }

        /*we now have the ke matrix*/
        addKe(e, ke);
    
        /*force vector*/
        double fe[4];

        for (int a=0;a<4;a++)
        {    
            /*second term int(na*h), always zero since support only h=0*/
            double fh=0;
            
            /*third term, -sum(kab*qb)*/
            double fg = 0;
            for (int b=0;b<4;b++)
            {
                int n = tet.con[b];
                double gb = boundary_potential[n];
                fg-=ke[a][b]*gb;
            }

            /*combine*/
            fe[a] = fh + fg;
        }

        addFe(F0, e,fe);
    }  /*end of cell*/


}

//*************************************************************************************************
/*computes "ff" part of F*/
void FESolver::buildF1Vector(const double *ion_den)
{
    double *f = new double[neq];
    /*start by computing the RHS term on all unknown nodes*/
    for (int n=0; n<n_nodes; n++)
    {
        int A = ID[n];
        if (A<0) continue;    /*skip known nodes*/
        f[A] = (QE/EPS0)*(ion_den[n]+n0*exp((p_d[A]-phi0)/kTe));
    }

    /*loop over cells*/
    for (int e=0; e<n_cells; e++)
    {
        double fe[4];
        for (int a=0; a<4; a++)
        {    
            /*first term is int(na*f), set to zero for now*/
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/
            {
                /*perform quadrature*/
                for (int k=0; k<n_int; k++)
                {
                    for (int j=0; j<n_int; j++)
                    {
                        for (int i=0; i<n_int; i++)
                        {
                            /*change of limits*/
                            double xi = 0.5*(p_l[i]+1);
                            double eta = 0.5*(p_l[j]+1);
                            double zeta = 0.5*(p_l[k]+1);

                            double Na = evalNa(a,xi,eta,zeta);

                            ff += f[A]*Na*detJ[e]*p_W[i]*p_W[j]*p_W[k];
                        }
                    }
                }
                ff*=(1.0/8.0);    /*change of limits*/
                fe[a] = ff;
            }
        }

        addFe(F1, e,fe);
    }

    delete[] f;
}

//*************************************************************************************************
/*simple Gauss-Seidel solver for A*x=b*/
void FESolver::solveLinear_petsc(double **p_A, double *p_x, double *p_b)
{ TRACE_ME; 

    // printf("At solveLinear_petsc neq %d\n", neq);

    auto t1 = std::chrono::system_clock::now();

    initialze_matrix(p_A);
    
    auto t2 = std::chrono::system_clock::now();
    
    VecSetValues(b, neq, vecCol, p_b, INSERT_VALUES);
    VecAssemblyBegin(b); VecAssemblyEnd(b);
    
    auto t3 = std::chrono::system_clock::now();
    
    KSPSolve(ksp, b, x);
    
    auto t4 = std::chrono::system_clock::now();
    
    VecGetValues(x, neq, vecCol, p_x);
    
    auto t5 = std::chrono::system_clock::now();
    
    // PetscPrintf(PETSC_COMM_WORLD, "PETSC : Norm of error %g iterations %" PetscInt_FMT "\n", (double)norm, its);

    std::chrono::duration<double> d1 = t2-t1;
    std::chrono::duration<double> d2 = t3-t2;
    std::chrono::duration<double> d3 = t4-t3;
    std::chrono::duration<double> d4 = t5-t4;
    std::chrono::duration<double> dt = t5-t1;
    std::cout << "solveLinear_petsc - Time t1: " << d1.count() << " t2: " << d2.count() << " t3: " << d3.count() << " t4: " << d4.count() <<  std::endl;
    std::cout << "solveLinear_petsc - Time <chrono>: " << dt.count() << " s" << std::endl;
}

//*************************************************************************************************
void FESolver::initialze_matrix(double **p_A)
{
    if (!matIndexCreated)
    {
        int max_fields = 0;
        for (int i=0; i<neq; i++)
        {
            std::vector<int> tempVec;    
            for (int j=0; j<neq; j++)
            {                
                if ((std::abs(p_A[i][j]) > 1e-12) || (i == j)) tempVec.push_back(j);            
            }

            matCol[i] = tempVec.size();
            matIndex[i] = new int[tempVec.size()];
            matIndexValues[i] = new double[tempVec.size()];

            std::copy(tempVec.begin(), tempVec.end(), matIndex[i]);

            max_fields = (max_fields > tempVec.size()) ? max_fields : tempVec.size();
        }

        MatSeqAIJSetPreallocation(A, (max_fields + 1), nullptr); // +1 just to be safe
        matIndexCreated = true;
    }

    for (int i=0; i<neq; i++)
    {
        for (int j=0; j<matCol[i]; j++) matIndexValues[i][j] = p_A[i][matIndex[i][j]];

        MatSetValues(A, 1, &i, matCol[i], matIndex[i], matIndexValues[i], INSERT_VALUES); 
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
}


//*************************************************************************************************
/*simple Gauss-Seidel solver for A*x=b*/
void FESolver::solveLinear(double **A, double *p_x, double *p_b)
{ TRACE_ME;
    int it;
    const double tol=1e-4;
    double L2;

    for (int u=0; u<neq; u++)
    {
        if (fabs(A[u][u]) < 1e-12)
        {
            cerr<<"FESolver : Zero diagonal on "<<u<<endl;
            // exit(-11);
        }        
    }

    bool converged=false;
    for (it=0; it<10000; it++)
    {
        for (int u=0; u<neq; u++)
        {
            /*skip over unused nodes*/
            if (fabs(A[u][u])<1e-12) continue;

            double sum=0;
            for (int v=0;v<neq;v++)
            {
                if (u==v) continue;
                sum+=A[u][v]*p_x[v];
            }
            p_x[u] = (p_b[u]-sum)/A[u][u];
        }

        /*periodically compute residue*/
        if (it%25 == 0)
        {
            double L=0;
            for (int u=0; u<neq; u++)
            {
                double sum=0;
                for (int v=0; v<neq; v++)
                    sum += A[u][v]*p_x[v];

                double r=p_b[u]-sum;
                L+=r*r;
            }

            L2 = sqrt(L)/neq;
            if (L2<tol) 
            {
                converged=true; 
                break;
            }
        }
    }

    if (!converged) 
    {
        cerr << "FESolver : GS failed to converge in " << it << " iterations, " << setprecision(3) << ": L2=" << L2 << endl;
        exit(-8);
    }
}

//*************************************************************************************************
void FESolver::SolveFields(
    oppic_dat ion_den_dat,
    oppic_dat field_potential_dat,
    oppic_dat boundary_potential_dat,
    oppic_dat electric_field_dat)
{
    const double *ion_den      = (double *)ion_den_dat->data;
    double *field_potential    = (double *)field_potential_dat->data;
    double *boundary_potential = (double *)boundary_potential_dat->data;
    double *electric_field     = (double *)electric_field_dat->data;

    /*call potential solver*/
    computePhi(ion_den, field_potential, boundary_potential); 

    updateElementElectricField(field_potential, electric_field);
}

//*************************************************************************************************
/*wrapper for solving the non-linear Poisson's equation, writing to field_potential*/
void FESolver::computePhi(const double *ion_den, double *field_potential, const double *boundary_potential)
{ TRACE_ME;
    /*solve the system*/
    solveNonLinear(ion_den);                            // simPIC field_solve_poissons_equation

    /*combine d and g to phi, should be similar to sum_laplace in simPIC*/
    // CAN BE WRITTEN TO OP_PAR_LOOP
    for (int n=0; n<n_nodes; n++)                // simPIC op_par_loop__field_solve_sum_laplace
    {
        /*zero on non-g nodes*/
        field_potential[n] = boundary_potential[n];    

        /*is this a non-g node?*/
        int A=ID[n];
        if (A>=0)
            field_potential[n] += p_d[A];
    }
}

//*************************************************************************************************
/*updates electric field at cells/ cells*/
// CAN BE WRITTEN TO OP_PAR_LOOP
void FESolver::updateElementElectricField(const double *field_potential, double *electric_field)        //simPIC field_solve_get_potential_gradient and op_par_loop__weight_fields_to_particles
{ TRACE_ME;
    /*interpolate electric field*/
    for (int cellID=0; cellID<n_cells; cellID++)
    {
        Tetra &tet = volume.cells[cellID];
        for (int d=0; d<3; d++) 
            electric_field[cellID * DIMENSIONS + d] = 0.0;

        for (int nodeCon=0; nodeCon<NODES_PER_CELL; nodeCon++)
        {
            int nodeID = tet.con[nodeCon];
            double nx[3];
            getNax(nx,cellID,nodeCon,0.5,0.5,0.5); // last 3 variables not used, copies derivative to nx at a given cell node

            /*minus sign since negative gradient*/
            for (int d=0;d<3;d++) 
                electric_field[cellID * DIMENSIONS + d] -= nx[d]*field_potential[nodeID];

            // if (cellID == 9249)
            // {
            //     printf("field_potential %d\t%+2.30lE\t%+2.30lE\t%+2.30lE\t%+2.30lE\n", nodeID, field_potential[nodeID], nx[0],nx[1],nx[2]);
            //     printf("%d\t%+2.30lE\t%+2.30lE\t%+2.30lE\n\n", nodeID, electric_field[cellID * DIMENSIONS + 0], electric_field[cellID * DIMENSIONS + 1],electric_field[cellID * DIMENSIONS + 2]);
            // }
        }
    }
}
