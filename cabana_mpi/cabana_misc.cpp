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

#include "cabana_defs.h"
#include "cabana_part_distribution.h"

using namespace opp;

class DataPointers;
void init_mesh(const Deck& deck, std::shared_ptr<DataPointers> m);
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
std::shared_ptr<DataPointers> load_mesh(const Deck& deck) {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    if (OPP_rank == OPP_ROOT)      
        init_mesh(deck, g_m);

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
void init_mesh(const Deck& deck, std::shared_ptr<DataPointers> m) {

    const OPP_INT nx             = deck.nx;
    const OPP_INT ny             = deck.ny;
    const OPP_INT nz             = deck.nz;
    const OPP_REAL c_widths[DIM] = { deck.dx, deck.dy, deck.dz };

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
void init_particles(const Deck& deck, opp_dat part_index, opp_dat part_pos, opp_dat part_vel, opp_dat part_streak_mid,
                    opp_dat part_weight, opp_dat part_mesh_rel, opp_dat cell_pos_ll, opp_dat cell_cgid) 
{
    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles START");

    const OPP_INT npart_per_cell = deck.nppc;

    const int cell_count = cell_pos_ll->set->size;
    const int rank_npart = npart_per_cell * cell_count;
    int rank_part_start  = 0;

    if (rank_npart <= 0) {
        opp_printf("Setup", "Error No particles to add in rank %d", OPP_rank);
        return;
    }

#ifdef USE_MPI // canculate the starting particle index incase of MPI
    {
        std::vector<OPP_INT> temp(OPP_comm_size, 0);
        MPI_Allgather(&rank_npart, 1, MPI_INT, temp.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < OPP_rank; ++i) rank_part_start += temp[i];
    }  
#endif

    if (OP_DEBUG)
        opp_printf("Setup", "%d parts to add in rank %d [part_start=%d]", rank_npart, OPP_rank, rank_part_start);

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles oppic_increase_particle_count rank_npart=%d", rank_npart);

    // Host/Device space to store the particles.
    oppic_increase_particle_count(part_index->set, rank_npart);

    if (opp_params->get<OPP_STRING>("part_enrich") == "two_stream")
        enrich_particles_two_stream(deck, cell_count, (OPP_REAL*)part_pos->data, (OPP_REAL*)part_vel->data, 
                (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data, (OPP_INT*)part_index->data, 
                (OPP_REAL*)part_weight->data, (OPP_REAL*)cell_pos_ll->data, rank_part_start, (OPP_INT*)cell_cgid->data);
    else
        enrich_particles_random(deck, cell_count, (OPP_REAL*)part_pos->data, (OPP_REAL*)part_vel->data, 
            (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data, (OPP_INT*)part_index->data, 
            (OPP_REAL*)part_weight->data, (OPP_REAL*)cell_pos_ll->data, rank_part_start);

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Uploading Start");

    opp_upload_particle_set(part_index->set);

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles END");
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


