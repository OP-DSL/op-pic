#include "simpic.h"
#include <cstring>

double dtfactor;

double *pdata;
double *lhsbuf, *rhsbuf;
int npart, lpart, rpart;

// fields arrays
double *Earray;
double *phiarray;
double *narray;
double *nback;

double area;
double density;
double np2c;
double q, m, qm;
double qscale;
double epsilon;
double wp;
double El, Er;
double dx, dt, t;
double xl, xr, L, Ll;
int ntimesteps;

int lproc, rproc;
int nproc, rank;
int ver, subver;
int last;

int ng, nc, ngglobal;

int nl, nr;
int nlp, nrp;  // counters for sending particle buffers
int nlr, nrr;

double *tri_a, *tri_b, *tri_c, *gam; // TODO : Enrich properly

bool diagnosticsflag;
double lhsvoltage;

double tttimes[NTIMES];
char names[NTIMES][STR_SIZE];
double tstart, tend;  // global times
double wtres;
double nparttot;

double tt;  // universal time

// variables for field communication
double *frbuffer;
double *fsbuffer;

char diags[NDIAGNOSTICS][STR_SIZE];
double * thediags[NDIAGNOSTICS];

char procname[STR_SIZE];


inline double xcoord(int i)
{
  return xl + i*dx;
}

void gradient(double *grad, double *arr, int n, double scale)
{
  int i;

  n--;
  for(i=1; i < n; i++)
    {
      grad[i] = scale*(arr[i+1] - arr[i-1]);
    }
}

void parsecmdline(int argc, char **argv)
{
    int i;
    int nnppcc;   // number of particles per cell
    int nnccpppp; // number of cells per processor 
    int nnnttt;   // number of time steps

    double ddttff, lhsdefault;

    ddttff = lhsdefault = 0.;
    nnppcc = nnccpppp = nnnttt = 0;

    for(i=1; i < argc; i++)
    {
        #ifdef DEBUG
        //      fprintf(stderr, "arg %d : %s\n", i, argv[i]);
        #endif

        if(strcmp(argv[i], "-ppc") == 0)
        {
            i++;
            sscanf(argv[i], "%d", &nnppcc);
        }
        else if(strcmp(argv[i], "-ncpp") == 0)
        {
            i++; 
            sscanf(argv[i], "%d", &nnccpppp);
        }
        else if(strcmp(argv[i], "-nt") == 0)
        {
            i++;
            sscanf(argv[i], "%d", &nnnttt);
        }
        else if(strcmp(argv[i], "-dtfactor") == 0)
        {
            i++;
            sscanf(argv[i], "%lf", &ddttff);
        }
        else if(strcmp(argv[i], "-lhsv") == 0)
        {
            i++;
            sscanf(argv[i], "%lf", &lhsdefault);
        }
        else // default case
        {
            fprintf(stderr, "\nError:\n");
            fprintf(stderr, "Unrecognized argument \"%s\"\n", argv[i]);
            //Quit();
        }
    }

    if((nnppcc < 1) || (nnccpppp < 1) || (nnnttt < 1))
    {
        fprintf(stderr, "\nError, input arguments must be entered!\n");
        //Quit();
    }

    if(ddttff <= 0.)
    {
        fprintf(stderr, "\nError, dtfactor MUST BE positive!\n");
        //Quit();
    }

    // set simulation variables from input data
    ntimesteps = nnnttt;
    nc = nnccpppp;
    ng = nc + 1;
    npart = nc*nnppcc;
    dtfactor = ddttff;
    lhsvoltage = lhsdefault;

}

void init(void)
{
    nparttot = 0.; // not used
    density = 1.E13;
    epsilon = 8.85E-12;
    area = 1.;
    L = 1.;
    q = 1.602E-19;
    m = 9.1E-31;

    rank = 0;
    nproc = 1;

    xl = rank/double(nproc);
    xr = (rank + 1.)/double(nproc);
    xl *= L;
    xr *= L;
    
    Ll = xr - xl;
    dx = Ll / double(nc);

    nl = rank - 1;
    nr = rank + 1;
    last = nproc - 1;
    if(rank == 0)
    {
        nl = last;
    }
    if(rank == last)
    {
        nr = 0;
    }

    np2c = density*area*Ll/double(npart);

    // calculate time step from plasma frequency
    wp = density*q*q/(epsilon*m);
    wp = sqrt(wp);

    dt = dtfactor/wp;

    qscale = np2c/(area*dx);

    t = ntimesteps*dt;
    qm = q*dt/m;

    diagnosticsflag = true;
}

double rhsV(double t)
{
    return 0.;
}

double lhsV(double t)
{
    return lhsvoltage;
}

void setcoeffs(double scale)
{
    if(rank == 0)
    {
        tri_a[0] = 0.0;
        tri_b[0] = 1.0;
        tri_c[0] = 0.0;
    }
    else
    {
        tri_b[0] = -scale*(1.0 + dx/xl);
        tri_c[0] = scale;
            
    }
    
    if(rank == last)
    {
        tri_a[nc] = 0.0;
        tri_b[nc] = 1.0;
        tri_c[nc] = 0.0;
    }
    else
    {
        tri_a[nc] = scale;
        tri_b[nc] = -scale*(1.0 + dx/(L - xr));
    }

    for(int i=1; i < nc; i++)
    {
        tri_a[i] = scale;
        tri_b[i] = (-2.0 * scale);
        tri_c[i] = scale;
    }
}


//*************************************************************************************************
void init_particles(
    double *&particle_position_x_tmp, 
    double *&particle_velocity_x_tmp, 
    double *&particle_field_E_tmp, 
    int *&particle_cell_index_tmp)
{
    particle_position_x_tmp    = new double[npart];
    particle_velocity_x_tmp    = new double[npart];
    particle_field_E_tmp    = new double[npart];
    particle_cell_index_tmp = new int[npart];

    for (int i=0; i < npart; i++)
    {
        double position             = xl + ((i + 1)*(xr-xl))/double(npart + 1);
        particle_position_x_tmp[i]  = position;
        particle_velocity_x_tmp[i]  = 0.0;
        particle_cell_index_tmp[i]  = ((position - xl)/dx);
        particle_field_E_tmp[i]     = 0.0;
    }
}


//*************************************************************************************************
void init_fields(
    double *&node_field_E_tmp, 
    double *&node_field_J_tmp, 
    double *&node_field_P_tmp,
    double *&node_xlocal_tmp, 
    int *&node_index_tmp, 
    int *&cell_to_nodes_tmp)
{
    node_field_E_tmp     = new double[ng];
    node_field_J_tmp     = new double[ng];
    node_field_P_tmp     = new double[ng];
    node_xlocal_tmp     = new double[ng];
    node_index_tmp         = new int[ng];
    cell_to_nodes_tmp     = new int[nc * 2];
    double xlocal = xl;

    for (int i=0, j=0; i < nc; i++, j+=2)
    {
        cell_to_nodes_tmp[j]     = i;
        cell_to_nodes_tmp[j+1]   = (i+1);
    }

    for (int i=0; i < ng; i++, xlocal +=dx)
    {
        node_field_E_tmp[i] = 0.0;      // Earray (Electric field)
        node_field_J_tmp[i] = 0.0;      // narray (Current density)
        node_field_P_tmp[i] = 0.0;      // phiarray (Potential)
        node_xlocal_tmp[i]  = xlocal;   // Local position of the node
        node_index_tmp[i]   = i;
    }

    {    
        tri_a   = new double[ng];
        tri_b   = new double[ng];
        tri_c   = new double[ng];
        gam     = new double[ng];

        for(int i=0; i < ng; i++)
        {
            tri_a[i] = tri_b[i] = tri_c[i] = 0.0;
        }

        setcoeffs(-epsilon/(q*dx*dx));
    }
} 


//*************************************************************************************************
void printDataToFiles(
    int n_nodes,
    int n_particles,
    double* node0_field_E,
    double* node0_field_J,
    double* node0_field_P,
    double* particle0_position_x,
    double* particle0_velocity_x,
    int* particle_cell_index,
    const std::string suffix)
{
    std::string fileNames[NDIAGNOSTICS] = {"_density.dat", "_E.dat", "_phi.dat", "_vxx.dat"};

    FILE *fp[NDIAGNOSTICS];
    for(int i = 0; i < NDIAGNOSTICS; i++)
    {    
        std::string file_name = std::string("files/") + suffix + fileNames[i]; 
        fp[i] = fopen(file_name.c_str(), "w");
    }
    
    for (int i = 0; i < n_nodes; i++)
    {
        fprintf(fp[0], "%d,%+2.30lE\n", i, node0_field_J[i]);
        fprintf(fp[1], "%d,%+2.30lE\n", i, node0_field_E[i]);
        fprintf(fp[2], "%d,%+2.30lE\n", i, node0_field_P[i]);
    }

    for (int i = 0; i < n_particles; i++)
    {
        fprintf(fp[3], "%d,%d,%+2.30lE,%+2.30lE\n", i, particle_cell_index[i], particle0_position_x[i], particle0_velocity_x[i]);
    }

    for(int i = 0; i < NDIAGNOSTICS; i++)
    {
        fclose(fp[i]);
    }
}

//*************************************************************************************************
void seq_field_solve_poissons_equation(
    int set_size,
    double* node0_field_J,
    double* node0_field_P)
{
    // modify density array
    double nlold, nrold;
    if(rank == 0)
    {
        nlold = node0_field_J[0];
        node0_field_J[0] = 0.;
    }
    else
    {
        node0_field_J[0] *= 2;
    }
    if(rank == last)
    {
        nrold = node0_field_J[nc];
        node0_field_J[nc] = 0.;
    }
    else
    {
        node0_field_J[nc] *= 2;
    }

    int nstrt = 0;

    // Tridiagonal matrix of Poisson equation 𝜙𝑗+1−2𝜙𝑗+𝜙𝑗−1=𝑝𝑗 is solved with Gaussian elimination by each processor
    {
        int j;
        double bet = tri_b[nstrt];
        node0_field_P[nstrt] = node0_field_J[nstrt]/bet;

        for(j = nstrt + 1; j < set_size; j++) 
        {
            gam[j]              = tri_c[j-1]/bet;
            bet                 = tri_b[j] - tri_a[j]*gam[j];
            node0_field_P[j]    = (node0_field_J[j] - tri_a[j]*node0_field_P[j-1])/bet;
        }

        for(j = nc - 1; j >= nstrt; j--)
        {
            node0_field_P[j] -= (gam[j+1] * node0_field_P[j+1]);
        }
    }

    // restore density array
    if(rank == 0)
    {
        node0_field_J[0] = 2*nlold;
    }
    if(rank == last)
    {
        node0_field_J[nc] = 2*nrold;
    }

}

//*************************************************************************************************
void seq_field_solve_get_potential_gradient(
    int set_size,
    double* nodes_field_E,
    const double* nodes_field_P)
{
    double scale = -0.5/dx;

    for(int i = 0; i < (set_size - 1); i++)
    {
        nodes_field_E[i] = scale * (nodes_field_P[i+1] - nodes_field_P[i-1]);
    }

    nodes_field_E[0]  = -(nodes_field_P[1] - nodes_field_P[0])/dx;
    nodes_field_E[nc] = -(nodes_field_P[nc] - nodes_field_P[nc-1])/dx;
}