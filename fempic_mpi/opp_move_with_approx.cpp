
#include "opp_move_with_approx.h"

using namespace opp;


#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

#define ASSIGN_CENTROID_TO_DIM(K)                                   \
    if (coordinate.K + this->gridSpacing <= maxCoordinate.K) {      \
        centroid.K = coordinate.K + this->gridSpacing / 2;;         \
    }                                                               \
    else {                                                          \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5;        \
    }                                                               \


//*******************************************************************************
CellApproximator::CellApproximator(const opp_dat nodePos_dat, double gridSpacing)
    : gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing) { 
    
    // if (OP_DEBUG)
        opp_printf("CellApproximator", "Create");

    generateBoundingBox(nodePos_dat);
    calculateGridDimensions();
} 

//*******************************************************************************
CellApproximator::~CellApproximator() { 
    
    // if (OP_DEBUG)
        opp_printf("CellApproximator", "Destroy");
}

//*******************************************************************************
opp_move_status CellApproximator::getCellIndexKernel(const double *point_pos, int* current_cell_index,
    double* point_lc, const double *cell_volume, const double *cell_det, const int *cell_connectivity) { 

    bool inside = true;  
    double coefficient2 = ONE_OVER_SIX / (*cell_volume);

    for (int i=0; i<N_PER_C; i++) { /*loop over vertices*/
    
        point_lc[i] = coefficient2 * (
            cell_det[i * DET_FIELDS + 0] - 
            cell_det[i * DET_FIELDS + 1] * point_pos[0] + 
            cell_det[i * DET_FIELDS + 2] * point_pos[1] - 
            cell_det[i * DET_FIELDS + 3] * point_pos[2]);
        
        if (point_lc[i] < 0.0 || 
            point_lc[i] > 1.0)  
                inside = false;
    }    
    
    if (inside) {
        return OPP_MOVE_DONE;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = point_lc[0];
    
    for (int i=1; i<NEIGHB_C; i++) {
        if (point_lc[i] < min_lc) {
            min_lc = point_lc[i];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i] >= 0) { // is there a neighbor in this direction?
        (*current_cell_index) = cell_connectivity[min_i];
        return OPP_NEED_MOVE;
    }
    else {
        (*current_cell_index) = MAX_CELL_INDEX;
        return OPP_NEED_REMOVE;
    }
}

//*******************************************************************************
const std::array<opp_point, 2>& CellApproximator::generateBoundingBox(const opp_dat node_pos_dat) {

    constexpr int DIM = 3;

    double* node_pos = (double*)node_pos_dat->data;

    for (int i = 0; i < node_pos_dat->set->size; i++) {
        minCoordinate.x = std::min(node_pos[i * DIM + 0], minCoordinate.x);
        minCoordinate.y = std::min(node_pos[i * DIM + 1], minCoordinate.y);
        minCoordinate.z = std::min(node_pos[i * DIM + 2], minCoordinate.z);
        maxCoordinate.x = std::max(node_pos[i * DIM + 0], maxCoordinate.x);
        maxCoordinate.y = std::max(node_pos[i * DIM + 1], maxCoordinate.y);
        maxCoordinate.z = std::max(node_pos[i * DIM + 2], maxCoordinate.z);
    }

    // if (OP_DEBUG)
        opp_printf("CellApproximator", "getBoundingBox Min[%2.20lE %2.20lE %2.20lE] Max[%2.20lE %2.20lE %2.20lE]", 
            minCoordinate.x, minCoordinate.y, minCoordinate.z, maxCoordinate.x, maxCoordinate.y, maxCoordinate.z);

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    return this->boundingBox;
}

//*******************************************************************************
const std::array<opp_point, 2>& CellApproximator::getBoundingBox() { 
    
    return this->boundingBox;
}

//*******************************************************************************
opp_point CellApproximator::getCentroidOfBox(const opp_point& coordinate) { 

    opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);

    constexpr int DIM = 3; // Hardcoded for now
    switch (DIM) {
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
            std::cerr << "Error getCentroidOfBox: Dimension invalid " << DIM << std::endl;
    }

    return centroid;
}

//*******************************************************************************
const std::vector<int>& CellApproximator::generateStructMeshToCellIndexVec(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
    const opp_map cellConnectivity_map) { 
    
    // if (OP_DEBUG)
        opp_printf("CellApproximator", "generateStructMeshToCellIndexVec");

    this->structMeshToCellIndexMap.clear();
    this->structMeshToCellIndexMap.reserve(gridDimensions.x * gridDimensions.y * gridDimensions.z);

    double lc[N_PER_C];
    int cell_set_size = cellVolume_dat->set->size;

    printf("grid dimensions %d %d %d\n", gridDimensions.x, gridDimensions.y, gridDimensions.z);

int ax =0, ay=0, az = 0, call = 0;
double x = 0.0, y = 0.0, z = 0.0;

    for (int dz = 0; dz < gridDimensions.z; dz++) { ax =0; ay=0; az++;
        
        z = minCoordinate.z + dz * this->gridSpacing;
        
        for (int dy = 0; dy < gridDimensions.y; dy++) { ax=0; ay++;
            
            y = minCoordinate.y + dy * this->gridSpacing;
            
            for (int dx = 0; dx < gridDimensions.x; dx++) { ax++;
                
                x = minCoordinate.x + dx * this->gridSpacing;
            
                const opp_point centroid = getCentroidOfBox(opp_point(x, y ,z));

                // Find in which cell this centroid lies and, in which MPI rank (for MPI backend)
                int cellIndex = 0;
                opp_move_status m = OPP_NEED_MOVE;

                do {
                    m = getCellIndexKernel(
                        (const double*)&centroid, 
                        &cellIndex,
                        lc,
                        &((double*)cellVolume_dat->data)[cellIndex], 
                        &((double*)cellDet_dat->data)[cellIndex * cellDet_dat->dim], 
                        &((int*)cellConnectivity_map->map)[cellIndex * cellConnectivity_map->dim]);

                } while (m == OPP_NEED_MOVE && cellIndex < cell_set_size);

                if (m == OPP_NEED_REMOVE) {
                    
                    // Eventhough the centroid is out of the structured mesh, 
                    // check atleast one vertex of the structured mesh can be within an unstructured mesh cell
                    const double gs = gridSpacing;
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
                        m = OPP_NEED_MOVE;

                        do {
                            m = getCellIndexKernel(
                                (const double*)&point, 
                                &cellIndex,
                                lc,
                                &((double*)cellVolume_dat->data)[cellIndex], 
                                &((double*)cellDet_dat->data)[cellIndex * cellDet_dat->dim], 
                                &((int*)cellConnectivity_map->map)[cellIndex * cellConnectivity_map->dim]);

                        } while (m == OPP_NEED_MOVE);

                        if (m == OPP_MOVE_DONE) {
                            found = true;
                            break;
                        }
                    }    

                    if (!found) {
                        cellIndex = MAX_CELL_INDEX;
                    }
                }

                size_t index1 = getCellIndexMappingIndex(opp_point(x, y, z));
                size_t index = (size_t)(dx + dy * gridDimensions.x + dz * gridDimensions.x * gridDimensions.y);
                if (index != this->structMeshToCellIndexMap.size()) {
                    {
                        size_t xIndex = static_cast<size_t>((x - this->minCoordinate.x) * this->oneOverGridSpacing);
                        size_t yIndex = static_cast<size_t>((y - this->minCoordinate.y) * this->oneOverGridSpacing);
                        size_t zIndex = static_cast<size_t>((z - this->minCoordinate.z) * this->oneOverGridSpacing);

                            // return (size_t)(xIndex + yIndex * gridDimensions.x + zIndex * gridDimensions.x * gridDimensions.y);
                        printf("There is an issue X %2.20lE %2.20lE %2.20lE | %zu %zu %zu | %d %d %d | %zu %zu %zu\n",
                        x,y,z,xIndex,yIndex,zIndex,ax, ay, az, index1, index,this->structMeshToCellIndexMap.size());
                    }
                    // {
                    //     size_t xIndex = static_cast<size_t>(((size_t)std::floor(x * 1e16) - (size_t)std::floor(this->minCoordinate.x * 1e16) ) * (this->oneOverGridSpacing * 1e-16));
                    //     size_t yIndex = static_cast<size_t>(((size_t)std::floor(y * 1e16) - (size_t)std::floor(this->minCoordinate.y * 1e16) ) * (this->oneOverGridSpacing * 1e-16));
                    //     size_t zIndex = static_cast<size_t>(((size_t)std::floor(z * 1e16) - (size_t)std::floor(this->minCoordinate.z * 1e16) ) * (this->oneOverGridSpacing * 1e-16));

                    //         // return (size_t)(xIndex + yIndex * gridDimensions.x + zIndex * gridDimensions.x * gridDimensions.y);
                    //     printf("There is an issue Y %2.20lE %2.20lE %2.20lE | %zu %zu %zu | %d %d %d | %zu %zu\n",
                    //     x,y,z,xIndex,yIndex,zIndex,ax, ay, az, index,this->structMeshToCellIndexMap.size());
                    // }
                }
                call++;
                this->structMeshToCellIndexMap.push_back(cellIndex);                
            }
        }
    }

    printf("CellApproximator grid dimensions %d %d %d -- calls %d %zu\n", ax, ay, az, call, this->structMeshToCellIndexMap.size());

    // printf("CellApproximator %d\n", (int)this->structMeshToCellIndexMap.size());
    // for (size_t i = 0; i < this->structMeshToCellIndexMap.size(); i++) {
    //     printf("%d ", this->structMeshToCellIndexMap[i]);

    //     if (i % 400 == 0) printf("\n");
    // }   
    // printf("\n");

    return this->structMeshToCellIndexMap;
}

//*******************************************************************************
const std::vector<int>& CellApproximator::getStructMeshToCellIndexVec() { 
    
    return this->structMeshToCellIndexMap;
}

//*******************************************************************************
bool CellApproximator::isCoordinateInBoundingBox(const opp_point& point) { 

    if (this->boundingBox[0].x > point.x || this->boundingBox[1].x < point.x) 
        return false;
    else if (this->boundingBox[0].y > point.y || this->boundingBox[1].y < point.y) 
        return false;
    else if (this->boundingBox[0].z > point.z || this->boundingBox[1].z < point.z) 
        return false;
    
    return true;
}

//*******************************************************************************
int CellApproximator::findClosestCellIndex(const opp_point& targetPosition) { 

    if (!isCoordinateInBoundingBox(targetPosition)) {
        // std::cerr << "isCoordinateInBoundingBox pos:" << targetPosition.x << "," << targetPosition.y << "," 
        //     << targetPosition.z << " is not within the bounding box" << std::endl;
        return MAX_CELL_INDEX;
    }

    int targetXIndex = static_cast<int>((targetPosition.x - minCoordinate.x) * this->oneOverGridSpacing);
    int targetYIndex = static_cast<int>((targetPosition.y - minCoordinate.y) * this->oneOverGridSpacing);
    int targetZIndex = static_cast<int>((targetPosition.z - minCoordinate.z) * this->oneOverGridSpacing);

    int closestIndex = targetXIndex + targetYIndex * gridDimensions.x + targetZIndex * gridDimensions.x * gridDimensions.y;
    
    // Assume closestIndex is within structMeshToCellIndexMap.size()
    return this->structMeshToCellIndexMap[closestIndex];
}

//*******************************************************************************
// This is correct, but has double precision issues
size_t CellApproximator::getCellIndexMappingIndex(const opp_point& position) {

    size_t xIndex = static_cast<size_t>((position.x - this->minCoordinate.x) * this->oneOverGridSpacing);
    size_t yIndex = static_cast<size_t>((position.y - this->minCoordinate.y) * this->oneOverGridSpacing);
    size_t zIndex = static_cast<size_t>((position.z - this->minCoordinate.z) * this->oneOverGridSpacing);

    return (size_t)(xIndex + yIndex * gridDimensions.x + zIndex * gridDimensions.x * gridDimensions.y);
}

//*******************************************************************************
const std::vector<opp_point>& CellApproximator::generateStructCoordinateVec() { 

    this->coordinateVec.clear(); 
    this->coordinateVec.reserve(gridDimensions.x * gridDimensions.y * gridDimensions.z);

    for (double z = minCoordinate.z; z < maxCoordinate.z; z += gridSpacing) {
        for (double y = minCoordinate.y; y < maxCoordinate.y; y += gridSpacing) {
            for (double x = minCoordinate.x; x < maxCoordinate.x; x += gridSpacing) {
                this->coordinateVec.emplace_back(opp_point(x, y, z)); 
            } 
        }
    }

    return this->coordinateVec;
}

//*******************************************************************************
const std::vector<opp_point>& CellApproximator::getStructuredCoordinateVec() { 

    return this->coordinateVec;
}

//*******************************************************************************
const opp_ipoint& CellApproximator::calculateGridDimensions() { 

    // double decimal point issues cause this to deviate
    // gridDimensions.x = (int)std::ceil((maxCoordinate.x - minCoordinate.x) * oneOverGridSpacing);
    // gridDimensions.y = (int)std::ceil((maxCoordinate.y - minCoordinate.y) * oneOverGridSpacing);
    // gridDimensions.z = (int)std::ceil((maxCoordinate.z - minCoordinate.z) * oneOverGridSpacing);

    // Hack : Use this to get correct grid dimension
    int ax =0, ay=0, az = 0, call = 0;
    for (double z = minCoordinate.z; z < maxCoordinate.z; z += this->gridSpacing) { 
        ax =0; ay=0; az++;
        for (double y = minCoordinate.y; y < maxCoordinate.y; y += this->gridSpacing) { 
            ax=0; ay++;
            for (double x = minCoordinate.x; x < maxCoordinate.x; x += this->gridSpacing) { 
                ax++;
            }
        }
    }

    gridDimensions.x = ax;
    gridDimensions.y = ay;
    gridDimensions.z = az;

    return gridDimensions;
}

//*******************************************************************************
const opp_ipoint& CellApproximator::getGridDimensions() { 

    return gridDimensions;
}

//*******************************************************************************
const std::vector<unsigned int>& CellApproximator::getHopCountsVec() {

    return this->hopCountsVec;
}

//*******************************************************************************
void CellApproximator::countHopsFromVec(const std::string& name, const std::vector<size_t>& vec) { 
    
    int less_2 = 0, less_3 = 0, less_4 = 0, less_5 = 0, less_10 = 0, less_50 = 0, less_100 = 0, less_500 = 0, 
        less_1000 = 0, more = 0, max = 0;

    for (auto& a : vec) {
        if (a < 2) less_2++;
        else if (a < 3) less_3++;
        else if (a < 4) less_4++;
        else if (a < 5) less_5++;
        else if (a < 10) less_10++;
        else if (a < 50) less_50++;
        else if (a < 100) less_100++;
        else if (a < 500) less_500++;
        else if (a < 1000) less_1000++;
        else more++;

        max = (max < a) ? a : max;
    }

    // if (OP_DEBUG)
        opp_printf("CellApproximator", "countHopsFromVec %s MAX=%d | <2=%d <3=%d <4=%d <5=%d <10=%d <50=%d <100=%d <500=%d <1000=%d >=1000=%d",
            name.c_str(), max, less_2, less_3, less_4, less_5, less_10, less_50, less_100, less_500, less_1000, more);
}

// #define USE_OMP
#ifdef USE_OMP

#include <omp.h>

//*******************************************************************************
void CellApproximator::move(const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat, const opp_dat cellVolume_dat, 
    const opp_dat cellDet_dat, const opp_map cellConnectivity_map) { 
    
    opp_profiler->start("MoveApprox");

    const double* pos = (const double*)pos_dat->data;
    int* cellIndex = (int*)cellIndex_dat->data;
    double* lc = (double*)lc_dat->data;
    opp_set set = pos_dat->set;
    int set_size = set->size;
    
    // hopCountsVec.resize(size);
    // std::fill(hopCountsVec.begin(), hopCountsVec.end(), 0); 

    opp_init_particle_move(set, 0, nullptr);

    int nthreads = omp_get_max_threads();

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;
        
        int cellIdx = MAX_CELL_INDEX;

        for (size_t i = start; i < finish; i++) {   
            
            opp_point* point = (opp_point*)&(pos[i * 3]);
            cellIndex[i] = findClosestCellIndex(*point);

            if (cellIndex[i] == MAX_CELL_INDEX) {
                part_remove_count_per_thr[thr] += 1;
                // std::cout << "Error... " << i << " " << cellIndex[i] << " " << cellVolume_dat->set->size << std::endl;
                continue;
            }

            opp_move_var m;

            do {
                cellIdx = cellIndex[i];

                m.move_status = getCellIndexKernel(
                    &((const double*)pos)[i * 3], 
                    &((int*)cellIndex)[i],
                    &((double*)lc)[i * 4],
                    &((double*)cellVolume_dat->data)[cellIdx], 
                    &((double*)cellDet_dat->data)[cellIdx * cellDet_dat->dim], 
                    &((int*)cellConnectivity_map->map)[cellIdx * cellConnectivity_map->dim]);
                
                // hopCountsVec[i]++;

            } while (opp_part_check_status(m, cellIdx, set, i, thr, thr));
        }
    }

    opp_finalize_particle_move(set);

    opp_profiler->end("MoveApprox");
}

#else

//*******************************************************************************
void CellApproximator::move(const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat, const opp_dat cellVolume_dat, 
    const opp_dat cellDet_dat, const opp_map cellConnectivity_map) { 
    
    opp_profiler->start("MoveApprox");

    const double* pos = (const double*)pos_dat->data;
    int* cellIndex = (int*)cellIndex_dat->data;
    double* lc = (double*)lc_dat->data;
    opp_set set = pos_dat->set;
    int size = set->size;
    int cellIdx = MAX_CELL_INDEX;

    hopCountsVec.resize(size);
    std::fill(hopCountsVec.begin(), hopCountsVec.end(), 0); 

    opp_init_particle_move(set, 0, nullptr);

    for (int i = 0; i < size; i++) {   
        
        opp_point* point = (opp_point*)&(pos[i * 3]);
        cellIndex[i] = findClosestCellIndex(*point);

        if (cellIndex[i] == MAX_CELL_INDEX) {
            set->particle_remove_count++;
            // std::cout << "Error... " << i << " " << cellIndex[i] << " " << cellVolume_dat->set->size << std::endl;
            continue;
        }

        opp_move_var m;

        do {
            cellIdx = cellIndex[i];

            m.move_status = getCellIndexKernel(
                &((const double*)pos)[i * 3], 
                &((int*)cellIndex)[i],
                &((double*)lc)[i * 4],
                &((double*)cellVolume_dat->data)[cellIdx], 
                &((double*)cellDet_dat->data)[cellIdx * cellDet_dat->dim], 
                &((int*)cellConnectivity_map->map)[cellIdx * cellConnectivity_map->dim]);
            
            hopCountsVec[i]++;

        } while (opp_part_check_status(m, cellIdx, set, i, set->particle_remove_count));
    }

    opp_finalize_particle_move(set);

    opp_profiler->end("MoveApprox");
}

#endif