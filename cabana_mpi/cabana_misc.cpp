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
void init_particles(std::shared_ptr<DataPointers> m);

int nx = -1;
int ny = -1;
int nz = -1;
int ng = -1;

//*************************************************************************************************
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
            if (c_ghost) delete[] c_ghost;
            if (c_iter_adv_e) delete[] c_iter_adv_e; 
            if (c_iter_acc) delete[] c_iter_acc;
           
            if (acc_cell_map) delete[] acc_cell_map;
            if (interp_cell_map) delete[] interp_cell_map;
            if (adv_b_cell_map) delete[] adv_b_cell_map;
            if (adv_e_cell_map) delete[] adv_e_cell_map;
            if (move_cell_map) delete[] move_cell_map;
            if (cell_cell_map) delete[] cell_cell_map;

            if (p_pos) delete[] p_pos;
            if (p_vel) delete[] p_vel;
            if (p_st_mid) delete[] p_st_mid;
            if (p_weight) delete[] p_weight;
            if (p_cid) delete[] p_cid;
            if (p_index) delete[] p_index;
          
            c_index = nullptr;
            c_e = nullptr;
            c_b = nullptr;
            c_j = nullptr;
            c_acc = nullptr;
            c_interp = nullptr;
            c_ghost = nullptr;
            c_iter_adv_e = nullptr;
            c_iter_acc = nullptr;
            
            acc_cell_map = nullptr;
            interp_cell_map = nullptr;
            adv_b_cell_map = nullptr;
            adv_e_cell_map = nullptr;
            move_cell_map = nullptr;
            cell_cell_map = nullptr;

            p_pos    = nullptr;
            p_vel    = nullptr;
            p_st_mid = nullptr;
            p_weight = nullptr;
            p_cid    = nullptr;
            p_index  = nullptr;
        }

        inline void CreateMeshArrays() {

            this->c_index      = new OPP_INT[this->n_cells * ONE];
            this->c_e          = new OPP_REAL[this->n_cells * DIM];  
            this->c_b          = new OPP_REAL[this->n_cells * DIM];        
            this->c_j          = new OPP_REAL[this->n_cells * DIM];
            this->c_acc        = new OPP_REAL[this->n_cells * DIM * ACCUMULATOR_ARRAY_LENGTH]; 
            this->c_interp     = new OPP_REAL[this->n_cells * INTERP_LEN]; 
            this->c_ghost      = new OPP_INT[this->n_cells * ONE];
            this->c_iter_adv_e = new OPP_INT[this->n_cells * ONE];
            this->c_iter_acc   = new OPP_INT[this->n_cells * ONE];

            this->acc_cell_map    = new OPP_INT[this->n_cells * 2*DIM];
            this->interp_cell_map = new OPP_INT[this->n_cells * 2*DIM];
            this->adv_b_cell_map  = new OPP_INT[this->n_cells * DIM];
            this->adv_e_cell_map  = new OPP_INT[this->n_cells * DIM];
            this->move_cell_map   = new OPP_INT[this->n_cells * 2*DIM];
            this->cell_cell_map   = new OPP_INT[this->n_cells * NEIGHBOUR_CELLS];
        }

        inline void CreateParticleArrays() {

            this->p_pos    = new OPP_REAL[this->n_particles * DIM];
            this->p_vel    = new OPP_REAL[this->n_particles * DIM];
            this->p_st_mid = new OPP_REAL[this->n_particles * DIM];
            this->p_weight = new OPP_REAL[this->n_particles];
            this->p_cid    = new OPP_INT[this->n_particles];
            this->p_index  = new OPP_INT[this->n_particles];
        }

        int n_cells     = 0;
        int n_particles = 0;

        OPP_INT* c_index         = nullptr;
        OPP_REAL* c_e            = nullptr;
        OPP_REAL* c_b            = nullptr;
        OPP_REAL* c_j            = nullptr;
        OPP_REAL* c_acc          = nullptr;
        OPP_REAL* c_interp       = nullptr;
        OPP_INT* c_ghost         = nullptr;
        OPP_INT* c_iter_adv_e    = nullptr;
        OPP_INT* c_iter_acc      = nullptr;
        OPP_INT* acc_cell_map    = nullptr;
        OPP_INT* interp_cell_map = nullptr;
        OPP_INT* adv_b_cell_map  = nullptr;
        OPP_INT* adv_e_cell_map  = nullptr;
        OPP_INT* move_cell_map   = nullptr;
        OPP_INT* cell_cell_map   = nullptr;

        OPP_REAL* p_pos    = nullptr;
        OPP_REAL* p_vel    = nullptr;
        OPP_REAL* p_st_mid = nullptr;
        OPP_REAL* p_weight = nullptr;
        OPP_INT*  p_cid    = nullptr;
        OPP_INT*  p_index  = nullptr;
};

//*************************************************************************************************
std::shared_ptr<DataPointers> LoadData() {

    std::shared_ptr<DataPointers> m(new DataPointers());

    init_mesh(m);
    // init_particles(m);

    opp_printf("LoadData", "Global counts - Cells[%d]", m->n_cells);

    return m;
}

//*************************************************************************************************
void init_mesh(std::shared_ptr<DataPointers> m) {

    nx = opp_params->get<OPP_INT>("nx");
    ny = opp_params->get<OPP_INT>("ny");
    nz = opp_params->get<OPP_INT>("nz");
    ng = opp_params->get<OPP_INT>("ng");
    // m->n_cells   = opp_params->get<OPP_INT>("num_cells");
    m->n_cells   = (nx + 2*ng) * (ny + 2*ng) * (nz + 2*ng);

	opp_printf("Setup", "init_mesh n_cells=%d nx=%d ny=%d nz=%d ng=%d", m->n_cells, nx, ny, nz, ng);

    m->CreateMeshArrays();

    // Enrich data
	for (int n = 0; n < m->n_cells; n++) {
        for (int d = 0; d < DIM; d++) {
            m->c_e[n * DIM + d] = 0.0;
            m->c_b[n * DIM + d] = 0.0;
            m->c_j[n * DIM + d] = 0.0;

            for (int i = 0; i < ACCUMULATOR_ARRAY_LENGTH; i++)
                m->c_acc[(n * DIM + d) * ACCUMULATOR_ARRAY_LENGTH + i] = 0.0;
        }

        for (int i = 0; i < INTERP_LEN; i++)
            m->c_interp[n * INTERP_LEN + i] = 0.0;
		
        m->c_index[n]      = n;
		m->c_ghost[n]      = 1;
		m->c_iter_adv_e[n] = 1;
		m->c_iter_acc[n]   = 0;
	}

	// Enrich mappings
	// for (int i = 0; i < 6 * m->n_cells; i++) {
	// 	m->interp_cell_map[i] = -1;
	// 	m->acc_cell_map[i] 	  = -1;
	// 	m->move_cell_map[i]	  = -1;
	// }
	// for (int i = 0; i < 3 * m->n_cells; i++) {
	// 	m->adv_b_cell_map[i] = -1;
	// 	m->adv_e_cell_map[i] = -1;
	// }
    for (int i = 0; i < NEIGHBOUR_CELLS * m->n_cells; i++) {
        m->cell_cell_map[i] = -1;
    }

	// size_t x_offset =  1;                       // VOXEL(x+1,y,  z,   nx,ny,nz);
    // size_t y_offset = (1*(nx+ng*2));            // VOXEL(x,  y+1,z,   nx,ny,nz);
    // size_t z_offset = (1*(nx+ng*2)*(ny+ng*2));  // VOXEL(x,  y,  z+1, nx,ny,nz);

	for (int x = 0; x < nx+2*ng; x++) {
		for (int y = 0; y < ny+2*ng; y++) {
			for (int z = 0; z < nz+2*ng; z++) {
                
                const int i = VOXEL(x,y,z, nx,ny,nz,ng);

                VOXEL_MAP(x-1, y-1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yd_zd]);
                VOXEL_MAP(x-1, y-1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yd_z]); 
                VOXEL_MAP(x-1, y-1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yd_zu]);
                VOXEL_MAP(x-1, y  , z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_y_zd]); 
                VOXEL_MAP(x-1, y  , z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_y_z]);  
                VOXEL_MAP(x-1, y  , z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_y_zu]); 
                VOXEL_MAP(x-1, y+1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yu_zd]);
                VOXEL_MAP(x-1, y+1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yu_z]); 
                VOXEL_MAP(x-1, y+1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xd_yu_zu]);
                VOXEL_MAP(x  , y-1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yd_zd]); 
                VOXEL_MAP(x  , y-1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yd_z]);  
                VOXEL_MAP(x  , y-1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yd_zu]); 
                VOXEL_MAP(x  , y  , z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_y_zd]);  
                VOXEL_MAP(x  , y  , z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_y_zu]);  
                VOXEL_MAP(x  , y+1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yu_zd]); 
                VOXEL_MAP(x  , y+1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yu_z]);  
                VOXEL_MAP(x  , y+1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::x_yu_zu]); 
                VOXEL_MAP(x+1, y-1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yd_zd]);
                VOXEL_MAP(x+1, y-1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yd_z]); 
                VOXEL_MAP(x+1, y-1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yd_zu]);
                VOXEL_MAP(x+1, y  , z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_y_zd]); 
                VOXEL_MAP(x+1, y  , z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_y_z]);  
                VOXEL_MAP(x+1, y  , z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_y_zu]); 
                VOXEL_MAP(x+1, y+1, z-1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yu_zd]);
                VOXEL_MAP(x+1, y+1, z  , nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yu_z]); 
                VOXEL_MAP(x+1, y+1, z+1, nx, ny, nz, ng, m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::xu_yu_zu]);

				if ( !(x < ng || y < ng || z < ng || x >= nx+ng || y >= ny+ng || z >= nz+ng) ) {

					m->c_ghost[i] = false;
					
                    // m->interp_cell_map[i * 6 + 0] = (i + x_offset);
					// m->interp_cell_map[i * 6 + 1] = (i + y_offset);
					// m->interp_cell_map[i * 6 + 2] = (i + z_offset);
					// m->interp_cell_map[i * 6 + 3] = (i + x_offset + y_offset);
					// m->interp_cell_map[i * 6 + 4] = (i + y_offset + z_offset);
					// m->interp_cell_map[i * 6 + 5] = (i + x_offset + z_offset);

                    // m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_0] = (i + x_offset);
					// m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_1] = (i + y_offset);
					// m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_2] = (i + z_offset);
					// m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_3] = (i + x_offset + y_offset);
					// m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_4] = (i + y_offset + z_offset);
					// m->cell_cell_map[i * NEIGHBOUR_CELLS + CellMap::interp_5] = (i + x_offset + z_offset);		
				
					// m->adv_b_cell_map[i * 3 + 0] = VOXEL(x+1, y,   z,   nx, ny, nz, ng);
					// m->adv_b_cell_map[i * 3 + 1] = VOXEL(x,   y+1, z,   nx, ny, nz, ng);
					// m->adv_b_cell_map[i * 3 + 2] = VOXEL(x,   y,   z+1, nx, ny, nz, ng);

					// m->move_cell_map[i * 6 + 0] = get_neighbour_cell(i, FACE_X_MIN , nx, ny, nz, ng);
					// m->move_cell_map[i * 6 + 1] = get_neighbour_cell(i, FACE_X_PLUS, nx, ny, nz, ng);
					// m->move_cell_map[i * 6 + 2] = get_neighbour_cell(i, FACE_Y_MIN , nx, ny, nz, ng);
					// m->move_cell_map[i * 6 + 3] = get_neighbour_cell(i, FACE_Y_PLUS, nx, ny, nz, ng);
					// m->move_cell_map[i * 6 + 4] = get_neighbour_cell(i, FACE_Z_MIN , nx, ny, nz, ng);
					// m->move_cell_map[i * 6 + 5] = get_neighbour_cell(i, FACE_Z_PLUS, nx, ny, nz, ng);	
                }

                // std::string log = "";
                // for (int d = 0; d < NEIGHBOUR_CELLS; d++)
                //     log += std::to_string(m->cell_cell_map[i * NEIGHBOUR_CELLS + d]) + " ";

                // opp_printf("INIT", "Cell %d | %d [%s] ", i, m->c_ghost[i], log.c_str());

				if (!(x < ng || y < ng || z < ng || x >= nx+2*ng || y >= ny+2*ng || z >= nz+2*ng)) {
					
                    if  ((VOXEL(x-1, y,   z,   nx,ny,nz,ng) > -1) &&
						 (VOXEL(x,   y-1, z,   nx,ny,nz,ng) > -1) &&
						 (VOXEL(x,   y,   z-1, nx,ny,nz,ng) > -1) &&
						 (VOXEL(x-1, y-1, z,   nx,ny,nz,ng) > -1) &&
						 (VOXEL(x,   y-1, z-1, nx,ny,nz,ng) > -1) &&
						 (VOXEL(x-1, y,   z-1, nx,ny,nz,ng) > -1)) {
                                
							// m->acc_cell_map[i * 6 + 0] = VOXEL(x-1, y,   z,   nx,ny,nz,ng);
							// m->acc_cell_map[i * 6 + 1] = VOXEL(x,   y-1, z,   nx,ny,nz,ng);
							// m->acc_cell_map[i * 6 + 2] = VOXEL(x,   y,   z-1, nx,ny,nz,ng);
							// m->acc_cell_map[i * 6 + 3] = VOXEL(x-1, y-1, z,   nx,ny,nz,ng);
							// m->acc_cell_map[i * 6 + 4] = VOXEL(x,   y-1, z-1, nx,ny,nz,ng);
							// m->acc_cell_map[i * 6 + 5] = VOXEL(x-1, y,   z-1, nx,ny,nz,ng);

							m->c_iter_acc[i] = true;	
						}	
				}

				if (x < ng || y < ng || z < ng) {
					m->c_iter_adv_e[i] = false;
				}
				// else {
				// 	m->adv_e_cell_map[i * 3 + 0] = VOXEL(x-1, y,   z,   nx, ny, nz, ng);
				// 	m->adv_e_cell_map[i * 3 + 1] = VOXEL(x,   y-1, z,   nx, ny, nz, ng);
				// 	m->adv_e_cell_map[i * 3 + 2] = VOXEL(x,   y,   z-1, nx, ny, nz, ng);				
				// }
			}
		}
	}
}

//*************************************************************************************************
void init_particles(std::shared_ptr<DataPointers> m) {

    const std::string file_name = opp_params->get<OPP_STRING>("part_file_name");
    m->n_particles              = opp_params->get<OPP_INT>("num_particles");

    FILE *fp;
	int size;
	double datd[7];
	int dati[1];

	if ((fp = fopen(file_name.c_str(), "r")) == NULL) {
		opp_printf("Setup", "Error can't open file %s", file_name.c_str());
		opp_abort();
	}
	if (fscanf(fp, "%d\n", &size) != 1) {
		opp_printf("Setup", "Error reading from %s", file_name.c_str());
		opp_abort();
	}

	opp_printf("Setup", "init_particles num_particles=%d size_from_file=%d", m->n_particles, size);

    m->CreateParticleArrays();

	for (int n = 0; n < m->n_particles; n++) {

		if (fscanf(fp, "%lf %lf %lf %lf %lf %lf %d %lf \n", 
            &datd[0], &datd[1], &datd[2], &datd[3], &datd[4], &datd[5], &dati[0], &datd[6]) != 8) {

			opp_printf("Setup", "Error reading from new_grid.dat for x");
			opp_abort();
		}
		else {
			m->p_pos[n * DIM + Dim::x] = datd[0];
			m->p_pos[n * DIM + Dim::y] = datd[1];
			m->p_pos[n * DIM + Dim::z] = datd[2];
            m->p_vel[n * DIM + Dim::x] = datd[3];
			m->p_vel[n * DIM + Dim::y] = datd[4];
			m->p_vel[n * DIM + Dim::z] = datd[5];
			m->p_cid[n]                = dati[0];
			m->p_weight[n]             = datd[6];
            m->p_index[n]              = n;

            m->p_st_mid[n * DIM + Dim::x] = 0.0;
			m->p_st_mid[n * DIM + Dim::y] = 0.0;
			m->p_st_mid[n * DIM + Dim::z] = 0.0;
		}
	}
}

// //*************************************************************************************************
// inline void init_particles(std::shared_ptr<DataPointers> m, const opp_dat cell_gcid, const opp_dat cell_ghost)
// {
//     opp_printf("Main", "init_particles ****");

//     const OPP_INT nx       = opp_params->get<OPP_INT>("nx");
//     const OPP_INT ny       = opp_params->get<OPP_INT>("ny");
//     const OPP_INT nz       = opp_params->get<OPP_INT>("nz");
//     const OPP_INT ng       = opp_params->get<OPP_INT>("ng");
//     const OPP_INT nppc     = opp_params->get<OPP_INT>("nppc");
//     const OPP_REAL w0      = opp_params->get<OPP_INT>("w0");
//     const OPP_REAL v0      = opp_params->get<OPP_INT>("v0");
//     const OPP_INT cellSize = 1.0;
//     m->n_particles         = nppc * nx * ny * nz;

//     // Calculate particle spacing
//     const double dx = cellSize / nppc;
//     const double dy = cellSize / nppc;
//     const double dz = cellSize / nppc;

//     m->CreateParticleArrays();

//     // create global cell index vs local index mapping
//     std::map<int, int> global_vs_local_cid;
//     for (int i = 0; i < cell_gcid->set->size; i++)
//     {
//         global_vs_local_cid.insert(std::pair<int, int>(((int*)cell_gcid->data)[i], i));
//     }

//     // Generate particles
//     int particle_index = 0;
//     for (int k = 0; k < nz; ++k) 
//     {
//         for (int j = 0; j < ny; ++j) 
//         {
//             for (int i = 0; i < nx; ++i) 
//             {
//                 const int gcell_index = VOXEL(i, j, k, nx,ny,nz,ng);
//                 auto it = global_vs_local_cid.find(gcell_index);
//                 if (it == global_vs_local_cid.end()) {
//                     opp_printf("Error", "gidx %d does not have a valid local cid", gcell_index);
//                     continue;
//                 }
//                 const int local_cell_index = it->second;

//                 if (((int*)cell_ghost->data)[local_cell_index] == 1)
//                     continue;

//                 for (int p = 0; p < x; ++p) 
//                 {
//                     m->p_pos[particle_index * DIM + Dim::x] = i * cellSize + p * dx;
//                     m->p_pos[particle_index * DIM + Dim::y] = j * cellSize + p * dy;
//                     m->p_pos[particle_index * DIM + Dim::z] = k * cellSize + p * dz;
//                     m->p_vel[particle_index * DIM + Dim::x] = 0.0;
//                     m->p_vel[particle_index * DIM + Dim::y] = v0;
//                     m->p_vel[particle_index * DIM + Dim::z] = 0.0;
//                     m->p_cid[particle_index]                = local_cell_index;
//                     m->p_weight[particle_index]             = w0;
//                     m->p_index[particle_index]              = particle_index;

//                     m->p_st_mid[particle_index * DIM + Dim::x] = 0.0;
//                     m->p_st_mid[particle_index * DIM + Dim::y] = 0.0;
//                     m->p_st_mid[particle_index * DIM + Dim::z] = 0.0;

//                     particle_index++;
//                 }
//             }
//         }
//     }
// }


//*************************************************************************************************
inline void init_particles(opp_dat part_index, opp_dat part_pos, opp_dat part_vel, opp_dat part_streak_mid, 
    opp_dat part_weight, opp_dat part_mesh_rel, const opp_dat cell_gcid, const opp_dat cell_ghost)
{
    opp_printf("Main", "init_particles ****");

    const OPP_INT nx        = opp_params->get<OPP_INT>("nx");
    const OPP_INT ny        = opp_params->get<OPP_INT>("ny");
    const OPP_INT nz        = opp_params->get<OPP_INT>("nz");
    const OPP_INT ng        = opp_params->get<OPP_INT>("ng");
    const OPP_INT nppc      = opp_params->get<OPP_INT>("nppc");
    const OPP_REAL w0       = opp_params->get<OPP_REAL>("part_weight");
    const OPP_REAL v0       = opp_params->get<OPP_REAL>("init_vel");
    const OPP_REAL c_len_x  = opp_params->get<OPP_REAL>("cell_len_x");
    const OPP_REAL c_len_y  = opp_params->get<OPP_REAL>("cell_len_y");
    const OPP_REAL c_len_z  = opp_params->get<OPP_REAL>("cell_len_z");
    
    opp_set part_set        = part_index->set;
    opp_set cell_set        = cell_ghost->set;

    // Calculate particle spacing
    const OPP_REAL dx = c_len_x / nppc;
    const OPP_REAL dy = c_len_y / nppc;
    const OPP_REAL dz = c_len_z / nppc;
    
    std::vector<OPP_INT> non_ghost_cells;
    for (int i = 0; i < cell_set->size; i++) 
    {
        if (((int*)cell_ghost->data)[i] != 1) 
            non_ghost_cells.push_back(i);
    }

    oppic_increase_particle_count(part_set, nppc * non_ghost_cells.size());

    int part_idx = 0;
    for (int i = 0; i < (int)non_ghost_cells.size(); i++) 
    {
        const int local_cid = non_ghost_cells[i];

        OPP_INT ix, iy, iz;
        RANK_TO_INDEX(((int*)cell_gcid->data)[local_cid], ix, iy, iz, (nx+(2*ng)), (ny+(2*ng)));

        for (int p = 0; p < nppc; ++p) 
        {
            ((OPP_REAL*)part_pos->data)[part_idx * DIM + Dim::x] = ix * c_len_x + p * dx;
            ((OPP_REAL*)part_pos->data)[part_idx * DIM + Dim::y] = iy * c_len_y + p * dy;
            ((OPP_REAL*)part_pos->data)[part_idx * DIM + Dim::z] = iz * c_len_z + p * dz;
            ((OPP_REAL*)part_vel->data)[part_idx * DIM + Dim::x] = 0.0;
            ((OPP_REAL*)part_vel->data)[part_idx * DIM + Dim::y] = v0;
            ((OPP_REAL*)part_vel->data)[part_idx * DIM + Dim::z] = 0.0;
            ((OPP_INT*)part_mesh_rel->data)[part_idx]            = local_cid;
            ((OPP_REAL*)part_weight->data)[part_idx]             = w0;

            // TODO: make correct part_index using previous non ghost cells and nppc
            ((OPP_INT*)part_index->data)[part_idx]               = part_idx; 

            ((OPP_REAL*)part_streak_mid->data)[part_idx * DIM + Dim::x] = 0.0;
            ((OPP_REAL*)part_streak_mid->data)[part_idx * DIM + Dim::y] = 0.0;
            ((OPP_REAL*)part_streak_mid->data)[part_idx * DIM + Dim::z] = 0.0;

            part_idx++;
        }       
    }

    opp_printf("Main", "%d particles added at init_particles ****", part_idx);
}

//*************************************************************************************************
inline void DistributeDataOverRanks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_particles), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<DataPointers>();

    m->n_cells     = opp_get_uniform_local_size(g_m->n_cells);
    m->n_particles = opp_get_uniform_local_size(g_m->n_particles);
    
    m->CreateMeshArrays();
    m->CreateParticleArrays();

    opp_uniform_scatter_array(g_m->c_index     , m->c_index     , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_e         , m->c_e         , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_b         , m->c_b         , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_j         , m->c_j         , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_acc       , m->c_acc       , g_m->n_cells, m->n_cells, DIM * ACCUMULATOR_ARRAY_LENGTH); 
    opp_uniform_scatter_array(g_m->c_interp    , m->c_interp    , g_m->n_cells, m->n_cells, INTERP_LEN); 
    opp_uniform_scatter_array(g_m->c_ghost     , m->c_ghost     , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_iter_adv_e, m->c_iter_adv_e, g_m->n_cells, m->n_cells, ONE); 
    opp_uniform_scatter_array(g_m->c_iter_acc  , m->c_iter_acc  , g_m->n_cells, m->n_cells, ONE); 

    opp_uniform_scatter_array(g_m->acc_cell_map   , m->acc_cell_map   , g_m->n_cells , m->n_cells, 2*DIM); 
    opp_uniform_scatter_array(g_m->interp_cell_map, m->interp_cell_map, g_m->n_cells , m->n_cells, 2*DIM); 
    opp_uniform_scatter_array(g_m->adv_b_cell_map , m->adv_b_cell_map , g_m->n_cells , m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->adv_e_cell_map , m->adv_e_cell_map , g_m->n_cells , m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->move_cell_map  , m->move_cell_map  , g_m->n_cells , m->n_cells, 2*DIM); 
    opp_uniform_scatter_array(g_m->cell_cell_map  , m->cell_cell_map  , g_m->n_cells , m->n_cells, NEIGHBOUR_CELLS); 

    opp_uniform_scatter_array(g_m->p_pos   , m->p_pos   , g_m->n_particles, m->n_particles, DIM); 
    opp_uniform_scatter_array(g_m->p_vel   , m->p_vel   , g_m->n_particles, m->n_particles, DIM); 
    opp_uniform_scatter_array(g_m->p_st_mid, m->p_st_mid, g_m->n_particles, m->n_particles, DIM); 
    opp_uniform_scatter_array(g_m->p_weight, m->p_weight, g_m->n_particles, m->n_particles, ONE);       
    opp_uniform_scatter_array(g_m->p_cid   , m->p_cid   , g_m->n_particles, m->n_particles, ONE); 
    opp_uniform_scatter_array(g_m->p_index , m->p_index , g_m->n_particles, m->n_particles, ONE); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif
}

//*************************************************************************************************
inline void genColoursForBlockPartition(opp_dat cell_colors, opp_dat cell_ghost, opp_dat cell_gbl_index)
{
    opp_printf("Main", "genColoursForBlockPartition **** %d", cell_colors->set->size);

#ifdef USE_MPI
    const OPP_INT nx           = opp_params->get<OPP_INT>("nx");
    const OPP_INT ny           = opp_params->get<OPP_INT>("ny");
    const OPP_INT nz           = opp_params->get<OPP_INT>("nz");
    const OPP_INT ng           = opp_params->get<OPP_INT>("ng");
    std::string log            = "";

    const OPP_INT* p_gbl_index = (OPP_INT*)cell_gbl_index->data;
    const OPP_INT* p_ghost     = (OPP_INT*)cell_ghost->data;
    OPP_INT* p_colors          = (OPP_INT*)cell_colors->data;

    std::vector<int> local_counts(OPP_comm_size);
    int old_value = 0;
    for (int i = 0; i < OPP_comm_size; i++) 
    {
        local_counts[i] = opp_get_uniform_local_size(nx * ny * nz, i) + old_value;
        old_value = local_counts[i];

        if (OPP_rank == OPP_ROOT) log += std::to_string(local_counts[i]) + " ";
    }
    if (OPP_rank == OPP_ROOT) opp_printf("COUNTS", "-> %s", log.c_str());

// global cell index - known ==> can find ix,iy,iz since nx,ny,nz,ng is known (RANK_TO_INDEX)
// can find whether it is a ghost or not, field is there
// Can find the uniform local sizes with nx,ny,nz and communicator size
// With ix,iy,iz and nx,ny,nz,ng=0 in VOXEL, we can find rank by comparing with uniform local sizes.
// what about ghost cells? iterate again and the ones without an MPI rank, put the neighbours MPI rank, 
// SIMPLE :) == Now implement!

    for (int i = 0; i < cell_colors->set->size; i++)
    {
        if (p_ghost[i] == 1)
        {
            p_colors[i] = OPP_rank;                                 // Changing this will reduce the halos
            // opp_printf("GHOST", "%d continue gbl %d", i, p_gbl_index[i]);
            continue;
        }

        OPP_INT ix, iy, iz;
        RANK_TO_INDEX(p_gbl_index[i], ix, iy, iz, (nx+(2*ng)), (ny+(2*ng)));

        const int non_ghost_index = VOXEL(ix-ng, iy-ng, iz-ng, nx,ny,nz,0);

        auto it = std::upper_bound(local_counts.begin(), local_counts.end(), non_ghost_index);
        if (it != local_counts.end()) 
        {
            p_colors[i] = std::distance(local_counts.begin(), it);
            // opp_printf("ASSIGN", "%d gbl %d assign %d", i, p_gbl_index[i], p_colors[i]);
        }
        else 
        {
            opp_printf("Main", "Error non_ghost_index %d cell_gbl_index %d", 
                non_ghost_index, p_gbl_index[i]);
            p_colors[i] = OPP_rank;
        }
    }
#endif
}

//*************************************************************************************************
inline void updatePartsWithLocalCellIndices(opp_dat cell_gcid, opp_dat part_cid, opp_dat part_index)
{
    opp_printf("Main", "updatePartsWithLocalCellIndices ****");

    std::map<int, int> global_vs_local_cid;

    for (int i = 0; i < cell_gcid->set->size; i++)
    {
        global_vs_local_cid.insert(std::pair<int, int>(((int*)cell_gcid->data)[i], i));
    }

    for (int i = 0; i < part_cid->set->size; i++)
    {
        auto it = global_vs_local_cid.find(((int*)part_cid->data)[i]);
        if (it == global_vs_local_cid.end()) {
            opp_printf("Error", "part lidx %d gidx %d with gcid %d does not have a valid local cid", 
                i, ((int*)part_index->data)[i], ((int*)part_cid->data)[i]);
            continue;
        }

        ((int*)part_cid->data)[i] = it->second;
    }
}

//*************************************************************************************************
inline void update_ghost_fields(const OPP_REAL* from, OPP_REAL* to)
{
    to[Dim::x] = from[Dim::x];
    to[Dim::y] = from[Dim::y];
    to[Dim::z] = from[Dim::z];
}

// Some values will not be in the current MPI rank, and VOXEL will give global cell id and not local cell id
// TODO : ERROR here for MPI
//*************************************************************************************************
inline void serial_update_ghosts_B(opp_dat cell_dat)
{
    if (FP_DEBUG) 
        opp_printf("CABANA", "serial_update_ghosts_B %s", cell_dat->name);

#ifdef USE_MPI
    int nargs = 1;
    oppic_arg args[nargs];
    args[0] = opp_get_arg(cell_dat,    OP_RW);
    args[0].idx = 2; // HACK to forcefully make halos to download

    opp_mpi_halo_exchanges_grouped(cell_dat->set, nargs, args, Device_CPU);
    opp_mpi_halo_wait_all(nargs, args);
#endif

    int from = -1, to = -1;
    const int dim = cell_dat->dim;

    // _zy_boundary
    for (int z = 1; z < nz+1; z++)
    {
        for (int y = 1; y < ny+1; y++)
        {
            from = VOXEL(1   , y, z, nx, ny, nz, ng);
            to   = VOXEL(nx+1, y, z, nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );

            // Copy x from RHS -> LHS
            from = VOXEL(nx  , y, z, nx, ny, nz, ng);
            to   = VOXEL(0   , y, z, nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );
        }
    }

    // _xz_boundary
    for (int x = 0; x < nx+2; x++)
    {
        for (int z = 1; z < nz+1; z++)
        {
            from = VOXEL(x,    1, z, nx, ny, nz, ng);
            to   = VOXEL(x, ny+1, z, nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );

            from = VOXEL(x, ny  , z, nx, ny, nz, ng);
            to   = VOXEL(x, 0   , z, nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );
        }
    }

    // _yx_boundary
    for (int x = 0; x < nx+2; x++)
    {
        for (int y = 0; y < ny+2; y++)
        {
            from = VOXEL(x, y, 1   , nx, ny, nz, ng);
            to   = VOXEL(x, y, nz+1, nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );

            from = VOXEL(x, y, nz  , nx, ny, nz, ng);
            to   = VOXEL(x, y, 0   , nx, ny, nz, ng);

            update_ghost_fields(
                &((double*)cell_dat->data)[from * dim], 
                &((double*)cell_dat->data)[to * dim]
            );
        }
    }

#ifdef USE_MPI
    opp_set_dirtybit(nargs, args);
#endif
}

//*************************************************************************************************
inline void serial_update_ghosts(opp_dat cell_j)
{
    if (FP_DEBUG) 
        opp_printf("CABANA", "serial_update_ghosts %s", cell_j->name);

#ifdef USE_MPI
    int nargs = 1;
    oppic_arg args[nargs];
    args[0] = opp_get_arg(cell_j,    OP_RW);
    args[0].idx = 2; // HACK to forcefully make halos to download
    
    MPI_Barrier(MPI_COMM_WORLD);

    opp_mpi_halo_exchanges_grouped(cell_j->set, nargs, args, Device_CPU);
    opp_mpi_halo_wait_all(nargs, args);

    MPI_Barrier(MPI_COMM_WORLD);
#endif

    double* p_cell_j = (double*)cell_j->data;
    const int dim = cell_j->dim;

    // _x_boundary
    for (int x = 0; x < nx+2; x++)
    {
        for(int z = 1; z <= nz+1; z++){
            //y first
            const int from = VOXEL(x, ny+1, z, nx, ny, nz, ng);
            const int to   = VOXEL(x, 1   , z, nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::x] += p_cell_j[from * dim + Dim::x];
        }

        for(int y = 1; y <= ny+1; y++){
            //z next
            const int from = VOXEL(x, y, nz+1, nx, ny, nz, ng);
            const int to   = VOXEL(x, y, 1   , nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::x] += p_cell_j[from * dim + Dim::x];
        }
    }

    //_y_boundary
    for (int y = 1; y < ny+1; y++)
    {
        for (int x = 1; x <= nx+1; x++){
            //z first
            const int from = VOXEL(x   , y, nz+1, nx, ny, nz, ng);
            const int to   = VOXEL(x   , y, 1   , nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::y] += p_cell_j[from * dim + Dim::y];
        }

        for (int z = 1; z <= nz+1; z++){
            //x next
            const int from = VOXEL(nx+1, y, z   , nx, ny, nz, ng);
            const int to   = VOXEL(1   , y, z   , nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::y] += p_cell_j[from * dim + Dim::y];
        }
    }

    //_z_boundary
    for (int z = 1; z < nz+1; z++)
    {
        for (int y = 1; y <= ny+1; y++){
            //x first
            const int from = VOXEL(nx+1, y   , z, nx, ny, nz, ng);
            const int to   = VOXEL(1   , y   , z, nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::z] += p_cell_j[from * dim + Dim::z];
        }

        for (int x = 1; x <= nx+1; x++){
            //y next
            const int from = VOXEL(x   , ny+1, z, nx, ny, nz, ng);
            const int to   = VOXEL(x   , 1   , z, nx, ny, nz, ng);

            p_cell_j[to * dim + Dim::z] += p_cell_j[from * dim + Dim::z];
        }
    }

#ifdef USE_MPI
    opp_set_dirtybit(nargs, args);
#endif
}

