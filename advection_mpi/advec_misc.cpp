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

#include "advec_defs.h"
#include <random>

using namespace opp;

class DataPointers;
void init_mesh(std::shared_ptr<DataPointers> m);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

//*************************************************************************************************
class DataPointers // This is just a placeholder for initializing // No use in DSL
{
    public:
        DataPointers() {}
        virtual ~DataPointers()
        {
            DeleteValues();   
        }

        inline void DeleteValues() {

            if (c_index) delete[] c_index;
            if (c_pos_ll) delete[] c_pos_ll;
            if (c_colors) delete[] c_colors;
            if (cell_cell_map) delete[] cell_cell_map;

            c_index = nullptr;
            c_pos_ll = nullptr;
            c_colors = nullptr;
            cell_cell_map = nullptr;
        }

        inline void CreateMeshArrays() {

            this->c_index       = new OPP_INT[this->n_cells * ONE];
            this->c_pos_ll      = new OPP_REAL[this->n_cells * DIM];  
            this->c_colors      = new OPP_INT[this->n_cells * ONE];        
            this->cell_cell_map = new OPP_INT[this->n_cells * NEIGHBOURS];

            for (int i = 0; i < this->n_cells * NEIGHBOURS; i++)
                this->cell_cell_map[i] = -1;
        }

        int n_cells     = 0;

        OPP_INT* c_index       = nullptr;
        OPP_REAL* c_pos_ll   = nullptr;
        OPP_INT* c_colors      = nullptr;
        OPP_INT* cell_cell_map = nullptr;
};

//*************************************************************************************************
std::shared_ptr<DataPointers> LoadData() {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    if (OPP_rank == OPP_ROOT)      
        init_mesh(g_m);

    std::shared_ptr<DataPointers> m;
    distribute_data_over_ranks(g_m, m);
    
    return m;
}

//*************************************************************************************************
// TODO : Can move this to Core
inline std::vector<int> get_local_cell_count_array(const int num_cells)
{
    std::vector<int> local_counts(OPP_comm_size);

#ifdef USE_MPI
    int old_value = 0;
    for (int i = 0; i < OPP_comm_size; i++) {
        local_counts[i] = opp_get_uniform_local_size(num_cells, i) + old_value;
        old_value = local_counts[i];
    }
#else
    local_counts[0] = num_cells;
#endif

    return local_counts;
}

//*************************************************************************************************
void init_mesh(std::shared_ptr<DataPointers> m) {

    OPP_INT nx = opp_params->get<OPP_INT>("nx");
    OPP_INT ny = opp_params->get<OPP_INT>("ny");
    OPP_REAL cell_width = opp_params->get<OPP_REAL>("cell_width");

    m->n_cells = (nx * ny);
    const std::vector<int> local_counts = get_local_cell_count_array(m->n_cells);

    m->CreateMeshArrays();

	opp_printf("Setup", "init_mesh n_cells=%d nx=%d ny=%d", m->n_cells, nx, ny);

	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
 
            const int i = VOXEL(x,y, nx);
            
            m->c_index[i] = i;
            m->c_pos_ll[i*DIM + Dim::x] = x * cell_width;
            m->c_pos_ll[i*DIM + Dim::y] = y * cell_width;

            const auto it = std::upper_bound(local_counts.begin(), local_counts.end(), i);
            if (it != local_counts.end()) {
                m->c_colors[i] = std::distance(local_counts.begin(), it);
            }
            else {
                opp_printf("Main", "Error unable to color cell %d [%d,%d]", i, x, y);
                m->c_colors[i] = OPP_rank;
            } 

            VOXEL_MAP(x-1, y,   nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_y]);
            VOXEL_MAP(x+1, y,   nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_y]);
            VOXEL_MAP(x,   y-1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yd]);
            VOXEL_MAP(x,   y+1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yu]);

            // VOXEL_MAP(x-1, y-1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yd]);
            // VOXEL_MAP(x-1, y+1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_yu]);
            // VOXEL_MAP(x+1, y-1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yd]);
            // VOXEL_MAP(x+1, y+1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_yu]);
        }
    }

    // TODO : Sanity check : Ideally there should not be any -1 mappings
}

//*************************************************************************************************
template <typename RNG>
inline std::vector<std::vector<double>>
    get_normal_distribution(const int N, const int ndim, const double mu, const double sigma, RNG &rng) {

    std::normal_distribution<> d{mu, sigma};
    std::vector<std::vector<double>> array(ndim);
    for (int dimx = 0; dimx < ndim; dimx++) {
        array[dimx] = std::vector<double>(N);
        for (int px = 0; px < N; px++) {
            array[dimx][px] = d(rng);
        }
    }

    return array;
}

//*************************************************************************************************
template <typename RNG>
inline std::vector<std::vector<double>>
    get_uniform_within_extents(const int N, const int ndim, const double *extents, RNG &rng) {

    std::uniform_real_distribution<double> uniform_rng(0.0, 1.0);
    std::vector<std::vector<double>> positions(ndim);

    for (int dimx = 0; dimx < ndim; dimx++) {
        positions[dimx] = std::vector<double>(N);
        const double ex = extents[dimx];
        for (int px = 0; px < N; px++) {
            positions[dimx][px] = ex * uniform_rng(rng);
        }
    }

    return positions;
}

//*************************************************************************************************
inline void uniform_within_cartesian_cells(const OPP_REAL* cell_pos_ll, const OPP_INT cells_set_size,
    const OPP_INT npart_per_cell, std::vector<std::vector<double>> &positions, std::vector<int> &cells,
    std::mt19937 rng) {
 
    OPP_REAL cell_width = opp_params->get<OPP_REAL>("cell_width");
    double extents[DIM] = { cell_width, cell_width };

    const int npart_total = npart_per_cell * cells_set_size;
    
    cells.resize(npart_total);

    positions.resize(DIM);
    positions[0] = std::vector<double>(npart_total);
    positions[1] = std::vector<double>(npart_total);

    for (int cx = 0; cx < cells_set_size; cx++) {

        const int index_start = cx * npart_per_cell;
        const int index_end = (cx + 1) * npart_per_cell;

        auto positions_ref_cell = get_uniform_within_extents(npart_per_cell, DIM, extents, rng);

        int index = 0;
        for (int ex = index_start; ex < index_end; ex++) {

            cells.at(ex) = cx;

            for (int dx = 0; dx < DIM; dx++) {
                positions.at(dx).at(ex) =
                    cell_pos_ll[cx * DIM + dx] + positions_ref_cell.at(dx).at(index);
            }
            index++;
        }
    }
}

//*************************************************************************************************
void init_particles(opp_dat part_index, opp_dat part_pos, opp_dat part_vel, opp_dat part_mesh_rel, 
                    opp_dat cell_pos_ll) 
{
    OPP_INT nx          = opp_params->get<OPP_INT>("nx");
    OPP_INT ny          = opp_params->get<OPP_INT>("ny");
    OPP_INT n_particles = opp_params->get<OPP_INT>("n_particles");

    std::mt19937 rng_pos(52234234 + OPP_rank);
    std::mt19937 rng_vel(52234231 + OPP_rank);

    const int cell_count = cell_pos_ll->set->size;
    const int global_cell_count = nx * ny;
    const int npart_per_cell = std::round((double) n_particles / (double) global_cell_count);
    const int rank_npart = npart_per_cell * cell_count;

    if (rank_npart <= 0) {
        opp_printf("Setup", "Error No particles to add in rank %d", OPP_rank);
    }

    std::vector<std::vector<double>> positions;
    std::vector<int> cells;

    // Sample particles randomly in each local cell.
    uniform_within_cartesian_cells((OPP_REAL*)cell_pos_ll->data, cell_count, npart_per_cell, 
                                    positions, cells, rng_pos);

    // Sample some particle velocities.
    auto velocities = get_normal_distribution(rank_npart, DIM, 0.0, 0.5, rng_vel);

    // Host space to store the particles.
    oppic_increase_particle_count(part_index->set, rank_npart);

    // Populate the host space with particle data.
    for (int px = 0; px < rank_npart; px++) {
        
        ((OPP_REAL*)part_pos->data)[px * DIM + Dim::x] = positions.at(Dim::x).at(px);
        ((OPP_REAL*)part_pos->data)[px * DIM + Dim::y] = positions.at(Dim::y).at(px);
        ((OPP_REAL*)part_vel->data)[px * DIM + Dim::x] = velocities.at(Dim::x).at(px);
        ((OPP_REAL*)part_vel->data)[px * DIM + Dim::y] = velocities.at(Dim::y).at(px);

        ((OPP_INT*)part_mesh_rel->data)[px]            = cells.at(px);
        ((OPP_INT*)part_index->data)[px]               = px; 
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//*************************************************************************************************
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<DataPointers>();
    m->n_cells     = opp_get_uniform_local_size(g_m->n_cells);
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_index      , m->c_index      , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_pos_ll     , m->c_pos_ll     , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_colors     , m->c_colors     , g_m->n_cells, m->n_cells, ONE); 
    opp_uniform_scatter_array(g_m->cell_cell_map, m->cell_cell_map, g_m->n_cells, m->n_cells, NEIGHBOURS); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
    
#else
    m = g_m;
#endif
}
