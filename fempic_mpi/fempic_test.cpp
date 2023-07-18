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

#include "fempic.h"
#include <list>

constexpr double MAX_REAL = std::numeric_limits<double>::max();
constexpr double MIN_REAL = std::numeric_limits<double>::min();

constexpr int MAX_INT = std::numeric_limits<int>::max();
constexpr int MIN_INT = std::numeric_limits<int>::min();

constexpr double ONE_OVER_SIX = (1.0 / 6.0);

constexpr int TEST_ITER_COUNT = 100;

struct opp_point {
    opp_point(double _x, double _y, double _z) {
        x = _x; 
        y = _y;
        z = _z;
    };
    opp_point() { };

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

//*************************************************************************************************
opp_move_status find_point_in_mesh_cells(const double *point_pos, int* current_cell_index,
    const double *cell_volume, const double *cell_det, const int *cell_connectivity) {

    bool inside = true;  
    double lc[N_PER_C];
    double coefficient2 = ONE_OVER_SIX / (*cell_volume);

    for (int i=0; i<N_PER_C; i++) { /*loop over vertices*/
    
        lc[i] = coefficient2 * (
            cell_det[i * DET_FIELDS + 0] - 
            cell_det[i * DET_FIELDS + 1] * point_pos[0] + 
            cell_det[i * DET_FIELDS + 2] * point_pos[1] - 
            cell_det[i * DET_FIELDS + 3] * point_pos[2]);
        
        if (lc[i] < 0.0 || 
            lc[i] > 1.0)  
                inside = false;
    }    
    
    if (inside) {
        return OPP_MOVE_DONE;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = lc[0];
    
    for (int i=1; i<NEIGHB_C; i++) {
        if (lc[i] < min_lc) {
            min_lc = lc[i];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i] >= 0) { // is there a neighbor in this direction?
        (*current_cell_index) = cell_connectivity[min_i];
        return OPP_NEED_MOVE;
    }
    else {
        (*current_cell_index) = -999;
        return OPP_NEED_REMOVE;
    }
}

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
}

// For now, implement only 3D, MinCoordinate[x,y,z], MaxCoordinate[x,y,z]
std::array<opp_point, 2> getBoundingBox(opp_dat node_pos_dat)
{
    opp_point minCoordinate(MAX_REAL, MAX_REAL, MAX_REAL);
    opp_point maxCoordinate(MIN_REAL, MIN_REAL, MIN_REAL);

    double* node_pos = (double*)node_pos_dat->data;

    for (int i = 0; i < node_pos_dat->set->size; i++) {
        minCoordinate.x = std::min(node_pos[i * DIM + 0], minCoordinate.x);
        minCoordinate.y = std::min(node_pos[i * DIM + 1], minCoordinate.y);
        minCoordinate.z = std::min(node_pos[i * DIM + 2], minCoordinate.z);
        maxCoordinate.x = std::max(node_pos[i * DIM + 0], maxCoordinate.x);
        maxCoordinate.y = std::max(node_pos[i * DIM + 1], maxCoordinate.y);
        maxCoordinate.z = std::max(node_pos[i * DIM + 2], maxCoordinate.z);
    }

    opp_printf("getBoundingBox", "[%2.6lE %2.6lE %2.6lE] x [%2.6lE %2.6lE %2.6lE]", 
        minCoordinate.x, minCoordinate.y, minCoordinate.z, maxCoordinate.x, maxCoordinate.y, maxCoordinate.z);

    std::array<opp_point, 2> boundingBox = { minCoordinate, maxCoordinate };

    return boundingBox;
}

opp_point getCentroidOfBox(const opp_point& coordinate, const opp_point& maxCoordinate, double gridSpacing) {
    
#define ASSIGN_CENTROID_TO_DIM(K)       \
    if (coordinate.K + gridSpacing <= maxCoordinate.K) {  \
        centroid.K = coordinate.K + gridSpacing / 2;;            \
    }                                                               \
    else { \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5; \
    } \

    opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);
    
    // can remove this validation I guess! check further...
    if (coordinate.x >= maxCoordinate.x || coordinate.y >= maxCoordinate.y || coordinate.z >= maxCoordinate.z)
        return centroid;

    const int dim = 3; // Hardcoded for now
    switch (dim) {
        case 1:
            ASSIGN_CENTROID_TO_DIM(x); 
            break;
        case 2:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            break;
        case 3:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            ASSIGN_CENTROID_TO_DIM(z);
            break;
        default:
            std::cerr << "Error getCentroidOfBox: Dimension invalid " << dim << std::endl;
    }

#undef ASSIGN_CENTROID_TO_DIM
    
    return centroid;
}

#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

void generateStructMeshToCellIndexMap(const opp_point& minCoordinate, const opp_point& maxCoordinate, 
    double gridSpacing, const opp_point& gridDimensions, const opp_dat cell_volume_dat, const opp_dat cell_det_dat, 
    const opp_map cell_connectivity_map, std::vector<int>& structMeshToCellIndexMap) {
    
    structMeshToCellIndexMap.clear();
    structMeshToCellIndexMap.reserve(gridDimensions.x * gridDimensions.y * gridDimensions.z);

    for (double z = minCoordinate.z; z < maxCoordinate.z; z += gridSpacing) {
        for (double y = minCoordinate.y; y < maxCoordinate.y; y += gridSpacing) {
            for (double x = minCoordinate.x; x < maxCoordinate.x; x += gridSpacing) {
            
                const opp_point centroid = getCentroidOfBox(opp_point(x, y ,z), maxCoordinate, gridSpacing);

                // Find in which cell this centroid lies and in which MPI rank (for MPI backend)
                int cellIndex = 0;
                opp_move_status m = OPP_NEED_MOVE;

                do {
                    m = find_point_in_mesh_cells(
                        (const double*)&centroid, 
                        &cellIndex,
                        &((double*)cell_volume_dat->data)[cellIndex], 
                        &((double*)cell_det_dat->data)[cellIndex * cell_det_dat->dim], 
                        &((int*)cell_connectivity_map->map)[cellIndex * cell_connectivity_map->dim]);

                } while (m == OPP_NEED_MOVE);

                if (m == OPP_NEED_REMOVE) {
                    
                    // Eventhough the centroid is out of the structured mesh, 
                    // check atleast one vertex of the structured mesh is within an unstructured mesh cell
                    double gs = gridSpacing;
                    bool found = false;

                    std::array<opp_point,8> vertices = {
                        opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z)),
                        opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z)),
                        opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
                        opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z+gs)),
                        opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z)),
                        opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z)),
                        opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z+gs)),
                        opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
                    };

                    for (const auto& point : vertices) {
                        cellIndex = 0;
                        opp_move_status m = OPP_NEED_MOVE;

                        do {
                            m = find_point_in_mesh_cells(
                                (const double*)&point, 
                                &cellIndex,
                                &((double*)cell_volume_dat->data)[cellIndex], 
                                &((double*)cell_det_dat->data)[cellIndex * cell_det_dat->dim], 
                                &((int*)cell_connectivity_map->map)[cellIndex * cell_connectivity_map->dim]);

                        } while (m == OPP_NEED_MOVE);

                        if (m == OPP_MOVE_DONE) {
                            found = true;
                            break;
                        }
                    }    

                    if (!found) {
                        cellIndex = -999;
                    }
                }

                structMeshToCellIndexMap.push_back(cellIndex);
            }
        }
    }
}

#undef GET_VERT

void moveWithHops(size_t size, const double* dist, int* ci1, const opp_dat cell_volume_dat, const opp_dat cell_det_dat, 
                    const opp_map cell_connectivity_map) {

    for (size_t i = 0; i < size; i++) {    
        opp_move_status m = OPP_NEED_MOVE;

        do {
            int cellIndex = ci1[i];

            m = find_point_in_mesh_cells(
                &((const double*)dist)[i * 3], 
                &((int*)ci1)[i],
                &((double*)cell_volume_dat->data)[cellIndex], 
                &((double*)cell_det_dat->data)[cellIndex * cell_det_dat->dim], 
                &((int*)cell_connectivity_map->map)[cellIndex * cell_connectivity_map->dim]);

        } while (m == OPP_NEED_MOVE);
    }
}

int findClosestCellIndex(const double* targetPosition, const std::vector<int>& structMeshToCellIndexMap,  
                const opp_point& gridDimensions, double oneOverGridSpacing, const opp_point& minCoordinate) {

    int targetXIndex = static_cast<int>((targetPosition[0] - minCoordinate.x) * oneOverGridSpacing);
    int targetYIndex = static_cast<int>((targetPosition[1] - minCoordinate.y) * oneOverGridSpacing);
    int targetZIndex = static_cast<int>((targetPosition[2] - minCoordinate.z) * oneOverGridSpacing);

    int closestIndex = targetXIndex + targetYIndex * gridDimensions.x + targetZIndex * gridDimensions.x * gridDimensions.y;

    return structMeshToCellIndexMap[closestIndex];
}

void moveDirect(size_t size, const double* dist, int* ci2, const opp_dat cell_volume_dat, const opp_dat cell_det_dat, 
    const opp_map cell_connectivity_map, const std::vector<int>& structMeshToCellIndexMap, const opp_point& gridDimensions, 
    double gridSpacing, const opp_point& minCoordinate) {
    
    const double oneOverGridSpacing = (1.0 / gridSpacing);

    for (size_t i = 0; i < size; i++) {   

        ci2[i] = findClosestCellIndex(&((const double*)dist)[i * 3], structMeshToCellIndexMap, gridDimensions, 
                                        oneOverGridSpacing, minCoordinate);
        if (ci2[i] < 0 || ci2[i] >= cell_volume_dat->set->size) {
            ci2[i] = -999;
            continue;
        }

        opp_move_status m = OPP_NEED_MOVE;

        do {
            int cellIndex = ci2[i];

            m = find_point_in_mesh_cells(
                &((const double*)dist)[i * 3], 
                &((int*)ci2)[i],
                &((double*)cell_volume_dat->data)[cellIndex], 
                &((double*)cell_det_dat->data)[cellIndex * cell_det_dat->dim], 
                &((int*)cell_connectivity_map->map)[cellIndex * cell_connectivity_map->dim]);

        } while (m == OPP_NEED_MOVE);
    }
}


void generateCoordinateVec(const opp_point& minCoordinate, const opp_point& maxCoordinate, double gridSpacing,
                            std::vector<opp_point>& coordinateVec, opp_point& gridDimensions) {

    coordinateVec.clear(); 
    coordinateVec.reserve(gridDimensions.x * gridDimensions.y * gridDimensions.z);

    for (double z = minCoordinate.z; z < maxCoordinate.z; z += gridSpacing) {
        for (double y = minCoordinate.y; y < maxCoordinate.y; y += gridSpacing) {
            for (double x = minCoordinate.x; x < maxCoordinate.x; x += gridSpacing) {
                coordinateVec.emplace_back(opp_point(x, y, z)); 
            } 
        }
    }
}

void calculateGridDimensions(const opp_point& minCoordinate, const opp_point& maxCoordinate, double gridSpacing,
                                opp_point& gridDimensions) {

    const double oneOverGridSpacing = (1.0 / gridSpacing);

    gridDimensions.x = std::ceil((maxCoordinate.x - minCoordinate.x) * oneOverGridSpacing);
    gridDimensions.y = std::ceil((maxCoordinate.y - minCoordinate.y) * oneOverGridSpacing);
    gridDimensions.z = std::ceil((maxCoordinate.z - minCoordinate.z) * oneOverGridSpacing);
}

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << 
        "/ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_1.param" << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    {
        opp_profiler->start("Setup");

        std::string log     = "";

        std::shared_ptr<FieldPointers> g_m, m; // g_m - global mesh, m - local mesh
        g_m = std::make_shared<FieldPointers>();
        
        if (OPP_rank == OPP_ROOT)
        {
            // Load using the original FemPIC loaders and distribute
            g_m = LoadMesh();
            opp_printf("Main", "Global counts - Nodes[%d] Cells[%d] IFaces[%d]", 
                g_m->n_nodes, g_m->n_cells, g_m->n_ifaces);
        }

        DistributeMeshOverRanks(g_m, m);

        opp_set node_set         = opp_decl_mesh_set(m->n_nodes, "mesh_nodes");
        opp_set cell_set         = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_set iface_set        = opp_decl_mesh_set(m->n_ifaces, "inlet_faces_cells");
        opp_set particle_set     = opp_decl_part_set("particles", cell_set); 

        opp_map cell_v_nodes_map = opp_decl_mesh_map(cell_set,  node_set, N_PER_C,  m->c_to_n, "c_v_n_map");
        opp_map cell_v_cell_map  = opp_decl_mesh_map(cell_set,  cell_set, NEIGHB_C, m->c_to_c,  "c_v_c_map"); 
        opp_map iface_v_cell_map = opp_decl_mesh_map(iface_set, cell_set, ONE,      m->if_to_c, "if_v_c_map"); 

        opp_dat cell_det         = opp_decl_mesh_dat(cell_set, ALL_DET,     DT_REAL, m->c_det, "c_det");  
        opp_dat cell_volume      = opp_decl_mesh_dat(cell_set, ONE,         DT_REAL, m->c_vol, "c_volume");        
        opp_dat cell_ef          = opp_decl_mesh_dat(cell_set, DIM,         DT_REAL, m->c_ef,  "c_ef");
        opp_dat cell_shape_deriv = opp_decl_mesh_dat(cell_set, N_PER_C*DIM, DT_REAL, m->c_sd,  "c_shape_deri"); 
        // opp_dat cell_id          = opp_decl_mesh_dat(cell_set, ONE,         DT_INT,  m->c_id,  "c_id"); 
        opp_dat cell_colors      = opp_decl_mesh_dat(cell_set, ONE,         DT_INT,  m->c_col, "c_colors");

        opp_dat node_volume      = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_vol,     "n_vol");        
        opp_dat node_potential   = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_pot,     "n_potential");     
        opp_dat node_charge_den  = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_ion_den, "n_charge_den");
        opp_dat node_pos         = opp_decl_mesh_dat(node_set, DIM, DT_REAL, m->n_pos,     "n_pos");     
        opp_dat node_type        = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_type,    "n_type");
        opp_dat node_bnd_pot     = opp_decl_mesh_dat(node_set, ONE, DT_REAL, m->n_bnd_pot, "n_bnd_pot");
        // opp_dat node_id          = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_id,      "n_id"); 
        // opp_dat node_colors      = opp_decl_mesh_dat(node_set, ONE, DT_INT,  m->n_color,   "n_colors");

        opp_dat iface_v_norm  = opp_decl_mesh_dat(iface_set, DIM, DT_REAL, m->if_v_norm, "iface_v_norm");        
        opp_dat iface_u_norm  = opp_decl_mesh_dat(iface_set, DIM, DT_REAL, m->if_u_norm, "iface_u_norm"); 
        opp_dat iface_norm    = opp_decl_mesh_dat(iface_set, DIM, DT_REAL, m->if_norm,   "iface_norm");     
        opp_dat iface_area    = opp_decl_mesh_dat(iface_set, ONE, DT_REAL, m->if_area,   "iface_area");
        opp_dat iface_dist    = opp_decl_mesh_dat(iface_set, ONE, DT_INT,  m->if_dist,   "iface_dist");
        opp_dat iface_n_pos   = opp_decl_mesh_dat(iface_set, DIM, DT_REAL, m->if_n_pos,  "iface_n_pos"); 
        // opp_dat iface_id      = opp_decl_mesh_dat(iface_set, ONE, DT_INT,  m->if_id,     "iface_id"); 

        opp_dat part_position = opp_decl_part_dat(particle_set, DIM,     DT_REAL, nullptr, "part_position");
        opp_dat part_velocity = opp_decl_part_dat(particle_set, DIM,     DT_REAL, nullptr, "part_velocity");    
        opp_dat part_lc       = opp_decl_part_dat(particle_set, N_PER_C, DT_REAL, nullptr, "part_lc");
        opp_dat part_mesh_rel = opp_decl_part_dat(particle_set, ONE,     DT_INT,  nullptr, "part_mesh_rel", true);

        opp_set dummy_part_set   = opp_decl_part_set("dummy particles", cell_set); 
        opp_dat dummy_part_rand  = opp_decl_part_dat(dummy_part_set, 2, DT_REAL, nullptr, "dummy_part_rand");

        m->DeleteValues();
        
        int n_parts_to_inject = InitializeInjectDistributions(iface_dist, iface_area, dummy_part_rand);

        opp_profiler->end("Setup");

        auto boundingBox = getBoundingBox(node_pos); // index 0 is min, index 1 is max
        double gridSpacing = 6e-5;  // 60um grid spacing will get all unstructured mesh cells captured in the structured mesh

        opp_point gridDimensions;
        calculateGridDimensions(boundingBox[0], boundingBox[1], gridSpacing, gridDimensions);

        // TODO : ideally, without coordinateVec, we can get the coordinate by calculating using minCoordinate
        // Lets keep this for now
        std::vector<opp_point> coordinateVec; 
        generateCoordinateVec(boundingBox[0], boundingBox[1], gridSpacing, coordinateVec, gridDimensions);
        std::cout << "coordinateVec size=" << coordinateVec.size() << std::endl;

        std::vector<int> structMeshToCellIndexMap;
        generateStructMeshToCellIndexMap(boundingBox[0], boundingBox[1], gridSpacing, gridDimensions, cell_volume,
                            cell_det, cell_v_cell_map, structMeshToCellIndexMap);

        std::cout << "structMeshToCellIndexMap size=" << structMeshToCellIndexMap.size() << std::endl;

        std::cout << "structMeshToCellIndexMap0 " << structMeshToCellIndexMap[0] << std::endl;
        std::cout << "structMeshToCellIndexMap1 " << structMeshToCellIndexMap[1] << std::endl;

        { // sanity check to test whether all unstructured mesh cells are captured over structured mesh
            int count = 0;
            for (int i = 0; i < cell_set->size; i++)
            {
                auto it = std::find(structMeshToCellIndexMap.begin(), structMeshToCellIndexMap.end(), i);
                if (it == structMeshToCellIndexMap.end()) {
                    std::cout << "Element " << i << " not found in the list." << std::endl;
                    count++;
                }            
            }
            std::cout << count << " cells does not directly relate with structured mesh" << std::endl;
        }

        // Trying to load random file to generate positions and simulate both scenarios, move with hops and moveDirectly
        {
            int total_size = -1, fsize = -1, fdim = -1;
            FILE *fp = NULL;
            std::string rand_file_path = opp_params->get<OPP_STRING>("rand_file");
            
            if ((fp = fopen(rand_file_path.c_str(), "r")) == NULL) {
                opp_printf("InitializeInjectDistributions", "Unable to open file %s\n", 
                    rand_file_path.c_str());
                opp_abort();
            }
            if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2) {
                opp_printf("InitializeInjectDistributions", "Error reading file data from %s\n", 
                    rand_file_path.c_str());
                opp_abort();
            }
            total_size = fsize * fdim;
            double* dist = new double[fsize * 3];
            for (int n = 0; n < fsize; n++) {
                if (fscanf(fp, " %lf %lf\n", &dist[n * 2 + 0], &dist[n * 2 + 1]) != 2) 
                {
                    opp_printf("InitializeInjectDistributions", "Error reading from %s at index %d\n", 
                        rand_file_path.c_str(), n);
                    opp_abort();
                }
                dist[n * 2 + 2] = (dist[n * 2 + 0] + dist[n * 2 + 1]) * (boundingBox[1].z - boundingBox[0].z) + boundingBox[0].z;
                dist[n * 2 + 1] = (dist[n * 2 + 1]) * (boundingBox[1].y - boundingBox[0].y) + boundingBox[0].y;
                dist[n * 2 + 0] = (dist[n * 2 + 0]) * (boundingBox[1].x - boundingBox[0].x) + boundingBox[0].x;
            }
            fclose(fp);  

            std::vector<int> ci1(fsize);
            std::vector<int> ci2(fsize);

            opp_printf("Main", "moveWithHops START *************");
            opp_profiler->start("moveWithHops");
            for (int x = 0; x < TEST_ITER_COUNT; x++) {
                     
                std::fill(ci1.begin(), ci1.end(), 0);
                moveWithHops(fsize, dist, ci1.data(), cell_volume, cell_det, cell_v_cell_map);
            }
            opp_profiler->end("moveWithHops");
            opp_printf("Main", "moveWithHops DONE *************");

            opp_printf("Main", "moveDirect START *************");
            opp_profiler->start("moveDirect");
            for (int x = 0; x < TEST_ITER_COUNT; x++) {
                     
                std::fill(ci2.begin(), ci2.end(), 0);
                moveDirect(fsize, dist, ci2.data(), cell_volume, cell_det, cell_v_cell_map, structMeshToCellIndexMap,
                    gridDimensions, gridSpacing, boundingBox[0]);
            }
            opp_profiler->end("moveDirect");
            opp_printf("Main", "moveDirect DONE *************");

            opp_printf("Main", "moveWithHops2 START *************");
            opp_profiler->start("moveWithHops2");
            for (int x = 0; x < TEST_ITER_COUNT; x++) {
                     
                std::fill(ci1.begin(), ci1.end(), 0);
                moveWithHops(fsize, dist, ci1.data(), cell_volume, cell_det, cell_v_cell_map);
            }
            opp_profiler->end("moveWithHops2");
            opp_printf("Main", "moveWithHops2 DONE *************");

            opp_printf("Main", "moveDirect2 START *************");
            opp_profiler->start("moveDirect2");
            for (int x = 0; x < TEST_ITER_COUNT; x++) {
                     
                std::fill(ci2.begin(), ci2.end(), 0);
                moveDirect(fsize, dist, ci2.data(), cell_volume, cell_det, cell_v_cell_map, structMeshToCellIndexMap,
                    gridDimensions, gridSpacing, boundingBox[0]);
            }
            opp_profiler->end("moveDirect2");
            opp_printf("Main", "moveDirect2 DONE *************");

            // test whether the cell indices calculated by both ways are equal or not!
            int wrong = 0;
            for (size_t k = 0; k < ci1.size(); k++) {
                if (ci1[k] != ci2[k]) {
                    std::cout << "Index:" << k << " ci1:" << ci1[k] << " ci2:" << ci2[k] << std::endl;
                    wrong++;
                }
            }
            std::cout << "Lines with issues:" << wrong << std::endl;
        }
    }

    opp_exit();

    if (OPP_rank == OPP_ROOT) 
        opp_printf("Main", "opp_exit DONE *************XXXX");

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(cell_v_nodes_map  , f.c_str(), "cell_v_nodes_map.dat");
// opp_print_dat_to_txtfile(node_charge_den, f.c_str(), "node_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(cell_shape_deriv, "cell_shape_deriv.dat");