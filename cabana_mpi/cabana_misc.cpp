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
// AUTO GENERATED CODE
//*********************************************

#include "cabana_defs.h"

using namespace opp;

class DataPointers;
void init_mesh(std::shared_ptr<DataPointers> m);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

//*************************************************************************************************
/**
 * @brief Utility class to temporarily hold the mesh data until it is loaded by OP-PIC
 */
class DataPointers
{
    public:
        DataPointers() {}
        virtual ~DataPointers()
        {
            DeleteValues();   
        }

        inline void DeleteValues() {

            if (c_index) delete[] c_index;
            if (c_e) delete[] c_e;
            if (c_b) delete[] c_b;
            if (c_j) delete[] c_j;
            if (c_acc) delete[] c_acc;
            if (c_interp) delete[] c_interp;
            if (c_pos_ll) delete[] c_pos_ll;

            if (cell_cell_map) delete[] cell_cell_map;

            c_index = nullptr;
            c_e = nullptr;
            c_b = nullptr;
            c_j = nullptr;
            c_acc = nullptr;
            c_interp = nullptr;
            c_pos_ll = nullptr;

            cell_cell_map = nullptr;
        }

        inline void CreateMeshArrays() {

            this->c_index  = new OPP_INT[this->n_cells * ONE];
            this->c_e      = new OPP_REAL[this->n_cells * DIM];  
            this->c_b      = new OPP_REAL[this->n_cells * DIM];        
            this->c_j      = new OPP_REAL[this->n_cells * DIM];
            this->c_acc    = new OPP_REAL[this->n_cells * ACC_LEN]; 
            this->c_interp = new OPP_REAL[this->n_cells * INTERP_LEN]; 
            this->c_pos_ll = new OPP_REAL[this->n_cells * DIM];

            this->cell_cell_map   = new OPP_INT[this->n_cells * NEIGHBOURS];
        }

        int n_cells     = 0;
        int n_particles = 0;

        OPP_INT* c_index   = nullptr;
        OPP_REAL* c_e      = nullptr;
        OPP_REAL* c_b      = nullptr;
        OPP_REAL* c_j      = nullptr;
        OPP_REAL* c_acc    = nullptr;
        OPP_REAL* c_interp = nullptr;
        OPP_REAL* c_pos_ll = nullptr;

        OPP_INT* cell_cell_map   = nullptr;
};

//*************************************************************************************************
/**
 * @brief Initialize the rank specific mesh data to a DataPointers utility class shared pointer
 * @return std::shared_ptr<DataPointers>
 */
std::shared_ptr<DataPointers> LoadData() {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    if (OPP_rank == OPP_ROOT)      
        init_mesh(g_m);

    std::shared_ptr<DataPointers> m;
    distribute_data_over_ranks(g_m, m);

    return m;
}

//*************************************************************************************************
/**
 * @brief Initializes the mesh using 2D (nx,ny) and cell_width values in the config file
 *          Expect this to run only on the ROOT MPI rank
 * @param m std::shared_ptr<DataPointers> loaded with mesh data
 * @return (void)
 */
void init_mesh(std::shared_ptr<DataPointers> m) {

    const OPP_INT nx             = opp_params->get<OPP_INT>("nx");
    const OPP_INT ny             = opp_params->get<OPP_INT>("ny");
    const OPP_INT nz             = opp_params->get<OPP_INT>("nz");
    const OPP_REAL c_widths[DIM] = { opp_params->get<OPP_REAL>("c_width_x"),
                                   opp_params->get<OPP_REAL>("c_width_y"),
                                   opp_params->get<OPP_REAL>("c_width_z") };
    m->n_cells   = (nx * ny * nz);

    opp_printf("Setup", "init_mesh global n_cells=%d nx=%d ny=%d nz=%d", m->n_cells, nx, ny, nz);

    m->CreateMeshArrays();

    for (int n = 0; n < m->n_cells; n++) {

        m->c_index[n]      = n;

        for (int d = 0; d < DIM; d++) {
            m->c_e[n * DIM + d] = 0.0;
            m->c_b[n * DIM + d] = 0.0;
            m->c_j[n * DIM + d] = 0.0;

            for (int i = 0; i < ACCUMULATOR_ARRAY_LENGTH; i++)
                m->c_acc[(n * DIM + d) * ACCUMULATOR_ARRAY_LENGTH + i] = 0.0;
        }

        for (int i = 0; i < INTERP_LEN; i++)
            m->c_interp[n * INTERP_LEN + i] = 0.0;
        
        for (int i = 0; i < NEIGHBOURS; i++)
            m->cell_cell_map[n * NEIGHBOURS + i] = -1;
    }

    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                
                const int i = VOXEL(x,y,z, nx,ny,nz);

                m->c_pos_ll[i*DIM + Dim::x] = x * c_widths[Dim::x];
                m->c_pos_ll[i*DIM + Dim::y] = y * c_widths[Dim::y];
                m->c_pos_ll[i*DIM + Dim::z] = z * c_widths[Dim::z];

                VOXEL_MAP(x-1, y-1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yd_zd]);
                VOXEL_MAP(x-1, y-1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yd_z]); 
                VOXEL_MAP(x-1, y-1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yd_zu]);
                VOXEL_MAP(x-1, y  , z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_y_zd]); 
                VOXEL_MAP(x-1, y  , z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_y_z]);  
                VOXEL_MAP(x-1, y  , z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_y_zu]); 
                VOXEL_MAP(x-1, y+1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yu_zd]);
                VOXEL_MAP(x-1, y+1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yu_z]); 
                VOXEL_MAP(x-1, y+1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yu_zu]);
                VOXEL_MAP(x  , y-1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yd_zd]); 
                VOXEL_MAP(x  , y-1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yd_z]);  
                VOXEL_MAP(x  , y-1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yd_zu]); 
                VOXEL_MAP(x  , y  , z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_y_zd]);  
                VOXEL_MAP(x  , y  , z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_y_zu]);  
                VOXEL_MAP(x  , y+1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yu_zd]); 
                VOXEL_MAP(x  , y+1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yu_z]);  
                VOXEL_MAP(x  , y+1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yu_zu]); 
                VOXEL_MAP(x+1, y-1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yd_zd]);
                VOXEL_MAP(x+1, y-1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yd_z]); 
                VOXEL_MAP(x+1, y-1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yd_zu]);
                VOXEL_MAP(x+1, y  , z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_y_zd]); 
                VOXEL_MAP(x+1, y  , z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_y_z]);  
                VOXEL_MAP(x+1, y  , z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_y_zu]); 
                VOXEL_MAP(x+1, y+1, z-1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yu_zd]);
                VOXEL_MAP(x+1, y+1, z  , nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yu_z]); 
                VOXEL_MAP(x+1, y+1, z+1, nx, ny, nz, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yu_zu]);
            }
        }
    }

    // TODO : Sanity Check : Check whether any mapping is still negative or not

    opp_printf("Setup", "init_mesh DONE");
}

//*************************************************************************************************
/**
 * @brief Initializes the particles in to the particle dats in the arguments, using the cell_pos_ll dat
 *          Expect this to run on every MPI rank
 * @param part_index - opp_dat : Particle index relative to rank. TODO: make this global
 * @param part_pos - opp_dat : Particle 3D position (x,y,z)
 * @param part_vel - opp_dat : Particle 3D velocity (x,y,z)
 * @param part_streak_mid - opp_dat : Particle 3D temporary position (x,y,z)
 * @param part_weight - opp_dat : Particle weight
 * @param part_mesh_rel - opp_dat : Particle belonging cell index 
 * @param cell_pos_ll - opp_dat : Lower left 2D position coordicate of the cell
 * @return (void)
 */
void init_particles(opp_dat part_index, opp_dat part_pos, opp_dat part_vel, opp_dat part_streak_mid,
                    opp_dat part_weight, opp_dat part_mesh_rel, opp_dat cell_pos_ll) 
{
    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles START");

    const OPP_INT npart_per_cell = opp_params->get<OPP_INT>("num_part_per_cell");
    const OPP_REAL weight        = opp_params->get<OPP_REAL>("part_weight");
    const OPP_REAL init_vel      = opp_params->get<OPP_REAL>("init_vel");
    const double extents[DIM]    = { opp_params->get<OPP_REAL>("c_width_x"),
                                     opp_params->get<OPP_REAL>("c_width_y"),
                                     opp_params->get<OPP_REAL>("c_width_z") };

    std::mt19937 rng_pos(52234234 + OPP_rank);

    const int cell_count        = cell_pos_ll->set->size;
    const int rank_npart        = npart_per_cell * cell_count;
    int rank_part_start = 0;

    if (rank_npart <= 0) {
        opp_printf("Setup", "Error No particles to add in rank %d", OPP_rank);
    }

#ifdef USE_MPI // canculate the starting particle index incase of MPI
    {
        std::vector<OPP_INT> temp(OPP_comm_size, 0);
        MPI_Allgather(&rank_npart, 1, MPI_INT, temp.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < OPP_rank; ++i) rank_part_start += temp[i];
    }  
#endif

    if (OP_DEBUG)
        opp_printf("Setup", "%d particles to add in rank %d [part start idx = %d]", 
                    rank_npart, OPP_rank, rank_part_start);

    std::mt19937 rng(52234234);
    std::uniform_real_distribution<double> uniform_rng(0.0, 1.0);

    std::vector<double> uniform_dist_vec(npart_per_cell * DIM);
    for (int px = 0; px < npart_per_cell; px++) {
        for (int dimx = 0; dimx < DIM; dimx++) {
            uniform_dist_vec[px * DIM + dimx] = extents[dimx] * uniform_rng(rng);
        }
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles oppic_increase_particle_count rank_npart=%d", rank_npart);

    // Host/Device space to store the particles.
    oppic_increase_particle_count(part_index->set, rank_npart);

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Load data to dats Start");

    // std::string log = "";
    // for (const double& num : uniform_dist_vec)
    //     log += str(num, "%.15f ");
    // opp_printf("W", "%s", log.c_str());

    // Populate the host space with particle data.
    int p_idx = 0;
    for (int cx = 0; cx < cell_count; cx++) {
        for (int px = 0; px < npart_per_cell; px++) {
            
            ((OPP_REAL*)part_pos->data)[p_idx * DIM + Dim::x] = uniform_dist_vec[px * DIM + Dim::x] + 
                                                            ((OPP_REAL*)cell_pos_ll->data)[cx * DIM + Dim::x];
            ((OPP_REAL*)part_pos->data)[p_idx * DIM + Dim::y] = uniform_dist_vec[px * DIM + Dim::y] + 
                                                            ((OPP_REAL*)cell_pos_ll->data)[cx * DIM + Dim::y];
            ((OPP_REAL*)part_pos->data)[p_idx * DIM + Dim::z] = uniform_dist_vec[px * DIM + Dim::z] + 
                                                            ((OPP_REAL*)cell_pos_ll->data)[cx * DIM + Dim::z];

            ((OPP_REAL*)part_vel->data)[p_idx * DIM + Dim::x] = 0.0;
            ((OPP_REAL*)part_vel->data)[p_idx * DIM + Dim::y] = init_vel;
            ((OPP_REAL*)part_vel->data)[p_idx * DIM + Dim::z] = 0.0;

            ((OPP_REAL*)part_streak_mid->data)[p_idx * DIM + Dim::x] = 0.0;
            ((OPP_REAL*)part_streak_mid->data)[p_idx * DIM + Dim::y] = 0.0;
            ((OPP_REAL*)part_streak_mid->data)[p_idx * DIM + Dim::z] = 0.0;

            ((OPP_INT*)part_mesh_rel->data)[p_idx] = cx;
            ((OPP_INT*)part_index->data)[p_idx]    = (cx * npart_per_cell + px + rank_part_start); // this might not exactly match with MPI versions
            ((OPP_REAL*)part_weight->data)[p_idx]  = weight;

            p_idx++;
        }
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Uploading Start");

    opp_upload_particle_set(part_index->set);

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles END");
}

inline std::vector<OPP_INT> cabana_get_cells_per_dim() {
    
    std::vector<OPP_INT> arr = { 
        opp_params->get<OPP_INT>("nx"),
        opp_params->get<OPP_INT>("ny"),
        opp_params->get<OPP_INT>("nz") 
    };

    if (OP_DEBUG)
        opp_printf("N", "%d %d %d", arr[0], arr[1], arr[2]);

    return arr;
}

inline std::array<OPP_REAL, DIM> cabana_get_c_widths() {

    std::array<OPP_REAL, DIM> arr = { 
        opp_params->get<OPP_REAL>("c_width_x"),
        opp_params->get<OPP_REAL>("c_width_y"),
        opp_params->get<OPP_REAL>("c_width_z") 
    };

    if (OP_DEBUG)
        opp_printf("c_widths", "%2.25lE %2.25lE %2.25lE", arr[0], arr[1], arr[2]);

    return arr;
}

inline std::array<OPP_REAL, DIM> cabana_get_cdt_d() {

    const OPP_REAL c = opp_params->get<OPP_REAL>("c");
    const OPP_REAL dt = opp_params->get<OPP_REAL>("dt");

    std::array<OPP_REAL, DIM> arr = {
        (c * dt / opp_params->get<OPP_REAL>("c_width_x")),
        (c * dt / opp_params->get<OPP_REAL>("c_width_y")),
        (c * dt / opp_params->get<OPP_REAL>("c_width_z"))
    };
    
    if (OP_DEBUG)
        opp_printf("cdt_d", "%2.25lE %2.25lE %2.25lE", arr[0], arr[1], arr[2]);

    return arr;
}

inline OPP_REAL cabana_get_qdt_2mc() {
    OPP_REAL a = (opp_params->get<OPP_REAL>("qsp") * opp_params->get<OPP_REAL>("dt") / 
                    (2 * opp_params->get<OPP_REAL>("me") * opp_params->get<OPP_REAL>("c")));
    
    if (OP_DEBUG)
        opp_printf("qdt_2mc", "%2.25lE", a);

    return a;
}

inline std::array<OPP_REAL, DIM> cabana_get_p() {

    const OPP_REAL c = opp_params->get<OPP_REAL>("c");
    const OPP_REAL dt = opp_params->get<OPP_REAL>("dt");
    const OPP_REAL frac = 1.0f;

    std::array<OPP_REAL, DIM> arr = {
        (opp_params->get<OPP_INT>("nx")>0) ? (frac * c * dt / opp_params->get<OPP_REAL>("c_width_x")) : 0,
        (opp_params->get<OPP_INT>("ny")>0) ? (frac * c * dt / opp_params->get<OPP_REAL>("c_width_y")) : 0,
        (opp_params->get<OPP_INT>("nz")>0) ? (frac * c * dt / opp_params->get<OPP_REAL>("c_width_z")) : 0
    };
    
    if (OP_DEBUG)
        opp_printf("p", "%2.25lE %2.25lE %2.25lE", arr[0], arr[1], arr[2]);

    return arr;
}

inline std::array<OPP_REAL, DIM> cabana_get_acc_coef() {

    const OPP_REAL dt = opp_params->get<OPP_REAL>("dt");
    const OPP_REAL dx = opp_params->get<OPP_REAL>("c_width_x");
    const OPP_REAL dy = opp_params->get<OPP_REAL>("c_width_y");
    const OPP_REAL dz = opp_params->get<OPP_REAL>("c_width_z");

    std::array<OPP_REAL, DIM> arr = {
        0.25 / (dy * dz * dt),
        0.25 / (dz * dx * dt),
        0.25 / (dx * dy * dt)
    };
    
    if (OP_DEBUG)
        opp_printf("acc_coef", "%2.25lE %2.25lE %2.25lE", arr[0], arr[1], arr[2]);

    return arr;
}

inline OPP_REAL cabana_get_dt_eps0() {
    
    OPP_REAL a = (opp_params->get<OPP_REAL>("dt") / opp_params->get<OPP_REAL>("eps"));

    if (OP_DEBUG)
        opp_printf("dt_eps0", "%2.25lE", a);

    return a;
}

//*************************************************************************************************
/**
 * @brief This block distributes temporary DataPointers from ROOT rank to other ranks
 * @param g_m - Global mesh of temporary shared pointer of DataPointers, Root Rank should have data
 * @param m - rank specific block partitioned mesh of temporary shared pointer of DataPointers
 * @return (void)
 */
inline void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_particles), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<DataPointers>();

    m->n_cells     = opp_get_uniform_local_size(g_m->n_cells);
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_index , m->c_index , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_e     , m->c_e     , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_b     , m->c_b     , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_j     , m->c_j     , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_acc   , m->c_acc   , g_m->n_cells, m->n_cells, DIM * ACCUMULATOR_ARRAY_LENGTH); 
    opp_uniform_scatter_array(g_m->c_interp, m->c_interp, g_m->n_cells, m->n_cells, INTERP_LEN); 
    opp_uniform_scatter_array(g_m->c_pos_ll, m->c_pos_ll, g_m->n_cells, m->n_cells, DIM);

    opp_uniform_scatter_array(g_m->cell_cell_map  , m->cell_cell_map  , g_m->n_cells , m->n_cells, NEIGHBOURS); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif
}


