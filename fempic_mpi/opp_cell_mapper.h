#pragma once

#include "opp_bounding_box.h"
#include "opp_particle_mover_kernel.h"

#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

#define ASSIGN_CENTROID_TO_DIM(K)                                   \
    if (coordinate.K + this->gridSpacing <= maxCoordinate.K) {      \
        centroid.K = coordinate.K + this->gridSpacing * 0.5;         \
    }                                                               \
    else {                                                          \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5;        \
    }                                                               \

namespace opp {

    class CellMapper {
    
    public:

        //*******************************************************************************
        CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, 
            const std::shared_ptr<Comm> comm = nullptr) 
            : boundingBox(boundingBox), gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing), 
                minGlbCoordinate(boundingBox->getGlobalMin()), comm(comm) {
            
            const opp_point& minGblCoordinate = boundingBox->getGlobalMin();
            const opp_point& maxGblCoordinate = boundingBox->getGlobalMax();

            // this->globalGridDims.x = std::ceil((maxGblCoordinate.x - minGblCoordinate.x) * oneOverGridSpacing);
            // this->globalGridDims.y = std::ceil((maxGblCoordinate.y - minGblCoordinate.y) * oneOverGridSpacing);
            // this->globalGridDims.z = std::ceil((maxGblCoordinate.z - minGblCoordinate.z) * oneOverGridSpacing);   

            // Hack : Use this to get correct grid dimension
            int ax = 0, ay = 0, az = 0;
            double x, y, z;
            for (z = minGblCoordinate.z; z < maxGblCoordinate.z; z += this->gridSpacing) { 
                az++; 
            }
            for (y = minGblCoordinate.y; y < maxGblCoordinate.y; y += this->gridSpacing) { 
                ay++; 
            }
            for (x = minGblCoordinate.x; x < maxGblCoordinate.x; x += this->gridSpacing) { 
                ax++; 
            }
            this->globalGridDims.x = ax + 1;
            this->globalGridDims.y = ay + 1;
            this->globalGridDims.z = az + 1; 

            if (OPP_rank == OPP_ROOT)
                opp_printf("CellMapper", "Global Grid Size - [%d %d %d] gridSpacing [%2.10lE]", 
                    this->globalGridDims.x, this->globalGridDims.y, this->globalGridDims.z, this->gridSpacing); 
            
            const opp_point& minLocalCoordinate = boundingBox->getLocalMin();
            const opp_point& maxLocalCoordinate = boundingBox->getLocalMax();

            // Find the local ranks grid start indices
            ax = 0; ay = 0; az = 0;
            for (double z = minGblCoordinate.z; (z < minLocalCoordinate.z); z += this->gridSpacing) { 
                az++; 
            }
            for (double y = minGblCoordinate.y; (y < minLocalCoordinate.y); y += this->gridSpacing) { 
                ay++; 
            }
            for (double x = minGblCoordinate.x; (x < minLocalCoordinate.x); x += this->gridSpacing) { 
                ax++; 
            }
            this->localGridStart.x = ax == 0 ? ax : (ax - 1);
            this->localGridStart.y = ay == 0 ? ay : (ay - 1);
            this->localGridStart.z = az == 0 ? az : (az - 1);         

            // Find the local ranks grid end indices
            ax = 0; ay = 0; az = 0;
            for (double z = minGblCoordinate.z; (z <= maxLocalCoordinate.z); z += this->gridSpacing) { 
                az++; 
            }
            for (double y = minGblCoordinate.y; (y <= maxLocalCoordinate.y); y += this->gridSpacing) { 
                ay++; 
            }
            for (double x = minGblCoordinate.x; (x <= maxLocalCoordinate.x); x += this->gridSpacing) { 
                ax++; 
            }
            this->localGridEnd.x = this->globalGridDims.x == ax ? ax : (ax + 1);
            this->localGridEnd.y = this->globalGridDims.y == ay ? ay : (ay + 1);
            this->localGridEnd.z = this->globalGridDims.z == az ? az : (az + 1);     

            //if (OP_DEBUG)
                opp_printf("CellMapper", "Local Grid - Min[%d %d %d] Max[%d %d %d]", 
                    this->localGridStart.x, this->localGridStart.y, this->localGridStart.z, 
                    this->localGridEnd.x, this->localGridEnd.y, this->localGridEnd.z); 
        }

        //*******************************************************************************
        ~CellMapper() { 

#ifdef ENABLE_MPI
            CHECK(MPI_Win_free(&this->win_structMeshToCellMapping))
            this->structMeshToCellMapping = nullptr;

            CHECK(MPI_Win_free(&this->win_structMeshToRankMapping))
            this->structMeshToRankMapping = nullptr;           
#else
            delete[] this->structMeshToCellMapping;
#endif
        };

        //*******************************************************************************
        inline opp_point getCentroidOfBox(const opp_point& coordinate) { 

            opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);
            const opp_point& maxCoordinate = boundingBox->getGlobalMax();

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
        // Returns the global cell index
        inline size_t findStructuredCellIndex(const opp_point& pos) { 

            if (!boundingBox->isCoordinateInGlobalBoundingBox(pos)) {
                // std::cerr << "isCoordinateInGlobalBoundingBox pos:" << pos.x << "," << pos.y << "," 
                //     << pos.z << " is not within the bounding box" << std::endl;
                return MAX_CELL_INDEX;
            }

            return getCellIndexMappingIndex(pos);
        }

        //*******************************************************************************
        // Returns the global cell index
        inline int findClosestGlobalCellIndex(const size_t& structCellIdx) { 
            
            if (OP_DEBUG) {
                if (structCellIdx >= globalGridSize) {
                    opp_printf("findClosestGlobalCellIndex", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                        structCellIdx, globalGridSize);
                    return MAX_CELL_INDEX;
                }
            }
            
            return this->structMeshToCellMapping[structCellIdx];
        }

        //*******************************************************************************
        // Returns the global cell index
        inline int findClosestLocalCellIndex(const size_t& structCellIdx) { 
            
            const int globalCellIndex = findClosestGlobalCellIndex(structCellIdx);

            //if (OPP_comm_size == 1)
                return globalCellIndex;
                
            //return getLocalCellIndexFromGlobal(globalCellIndex);
        }

        //*******************************************************************************
        // Returns the rank of the cell
        inline int findClosestCellRank(const size_t& structCellIdx) { 

            if (OP_DEBUG) {
                if (structCellIdx >= globalGridSize) {
                    opp_printf("findClosestCellRank", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                        structCellIdx, globalGridSize);
                    return MAX_CELL_INDEX;
                }
            }

            return this->structMeshToRankMapping[structCellIdx];
        }

// Be careful when implementing MPI
// 1: structMeshToCellMapping should contain global indices - DONE
// 2: it should be populated one copy per shared memory (need to handle overlappings)
// 3: once populated per shared memory, need to share the details with inter nodes (can use inter node comm root) - DONE, CHECK

        //*******************************************************************************
        inline void generateStructMeshToGlobalCellMappings(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map) { 

            // generateStructMeshToGlobalCellMappings_multiHop(cellVolume_dat, cellDet_dat, global_cell_id_dat, cellConnectivity_map);
            generateStructMeshToGlobalCellMappings_allSearch(cellVolume_dat, cellDet_dat, global_cell_id_dat, cellConnectivity_map);
        }

        //*******************************************************************************
        inline void generateStructMeshToGlobalCellMappings_multiHop(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map) { 

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif            
            // if (OP_DEBUG)
                opp_printf("CellMapper", "generateStructMeshToGlobalCellMappings_multiHop global grid dimensions %d %d %d",
                    globalGridDims.x, globalGridDims.y, globalGridDims.z);

            //int globalGridSize = (globalGridDims.x * globalGridDims.y * globalGridDims.z);
            opp_set set = cellVolume_dat->set;

            createStructMeshMappingArrays();

            double lc[N_PER_C];
            int cellSetSizeIncHalo = set->size + set->exec_size + set->nonexec_size;

            // const opp_point& minCoordinate = boundingBox->getLocalMin(); // required for GET_VERT
            const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT

            double x = 0.0, y = 0.0, z = 0.0;

            auto multihop_mover = [&](const opp_point& point, int& cellIndex, opp_move_status& m, int dx, int dy, int dz) {

                cellIndex = set->size / 2;

                do {
                    m = getCellIndexKernel(
                        (const double*)&point, 
                        &cellIndex,
                        lc,
                        &((double*)cellVolume_dat->data)[cellIndex], 
                        &((double*)cellDet_dat->data)[cellIndex * cellDet_dat->dim], 
                        &((int*)cellConnectivity_map->map)[cellIndex * cellConnectivity_map->dim]);

                } while (m == OPP_NEED_MOVE && cellIndex < cellSetSizeIncHalo);  
            };

            const opp_point& minGlbCoordinate = boundingBox->getGlobalMin();

            std::map<size_t, opp_point> removedCoordinates;

            if (cellSetSizeIncHalo > 0) {

                // Step 1 : Get the centroids of the structured mesh cells and try to relate them to unstructured mesh cell indices
                for (int dz = this->localGridStart.z; dz <= this->localGridEnd.z; dz++) {
                    
                    z = minGlbCoordinate.z + dz * this->gridSpacing;
                    
                    for (int dy = this->localGridStart.y; dy <= this->localGridEnd.y; dy++) {
                        
                        y = minGlbCoordinate.y + dy * this->gridSpacing;
                        
                        for (int dx = this->localGridStart.x; dx <= this->localGridEnd.x; dx++) {
                            
                            x = minGlbCoordinate.x + dx * this->gridSpacing;
                            
                            const opp_point centroid = getCentroidOfBox(opp_point(x, y ,z));

                            opp_move_status m = OPP_NEED_MOVE;
                            int cellIndex = 0;

                            // Find in which cell this centroid lies and, in which MPI rank (for MPI backend)
                            multihop_mover(centroid, cellIndex, m, dx,dy,dz);

                            size_t index = (size_t)(dx + dy * globalGridDims.x + dz * globalGridDims.x * globalGridDims.y);

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "Iterating index %zu | %d %d %d | %2.6lE %2.6lE %2.6lE", index, dx,dy,dz, x,y,z);  

                            if (m == OPP_NEED_REMOVE) {
                                removedCoordinates.insert(std::make_pair(index, opp_point(x, y ,z)));
    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "Marking index %zu | %d %d %d | %2.6lE %2.6lE %2.6lE", index, dx,dy,dz, x,y,z);  
                                continue;
                            }

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "Writing index %zu | %d %d %d | cellIndex %d gbl %d | %s | %d", index, dx,dy,dz, cellIndex, 
    //     ((int*)global_cell_id_dat->data)[cellIndex], (cellIndex < cellSetSizeIncHalo) ? "YES" : "NO", set->size);  

                            // Allow neighbours to write on-behalf of the current rank, to reduce issues
                            if (cellIndex < cellSetSizeIncHalo) {
                                this->structMeshToCellMapping[index] = ((int*)global_cell_id_dat->data)[cellIndex];
#ifdef ENABLE_MPI
                                if (cellIndex < set->size)
                                    this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank
                                else {
                                    
                                    std::map<int, opp_particle_comm_data>* rankVsCommData = nullptr;
                                    for (auto& x : opp_part_comm_neighbour_data) {

                                        if (x.first->cells_set == set) {
                                            rankVsCommData = &(x.second);
                                            break;
                                        }
                                    }

                                    auto it = rankVsCommData->find(cellIndex);
                                    if (it == rankVsCommData->end())
                                        opp_abort(std::string("Cant find comm data of local cell"));

                                    // get the correct neighbour rank
                                    this->structMeshToRankMapping[index] = it->second.cell_residing_rank; 
                                }           
#endif
                            }
                        } 
                    }
                }
            }

            // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            // structMeshToCellMapping was initialized to MAX_CELL_INDEX, hence reducing to MIN to get the correct cell index
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToRankMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter

            // for (auto it = removedCoordinates.begin(); it != removedCoordinates.end(); ) {

            //     size_t removedIndex = it->first;
                
            //     if (this->structMeshToRankMapping[removedIndex] != MAX_CELL_INDEX) {
            //         it = removedCoordinates.erase(it); // This structured index is already written by another rank
            //         opp_printf("CellMapper", "Already written %d to struct index %zu", this->structMeshToRankMapping[removedIndex], removedIndex);
            //     } else {
            //         ++it; 
            //     }
            // }
#endif           

            printStructuredMesh("After centroid mappings");

            if (cellSetSizeIncHalo > 0) {

                // Step 3 : Iterate over all the NEED_REMOVE points, try to check whether atleast one vertex of the structured mesh can be within 
                //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
                for (auto& p : removedCoordinates) {

                    size_t index = p.first;
                    double& x = p.second.x;
                    double& y = p.second.y;
                    double& z = p.second.z;

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "removedCoordinates Iterating index %zu | %2.6lE %2.6lE %2.6lE", index, x,y,z);  

                    if (this->structMeshToCellMapping[index] == MAX_CELL_INDEX) {
                        
                        // still no one has claimed that this cell belongs to it

                        const double gs = gridSpacing;
                        int mostSuitableCellIndex = MAX_CELL_INDEX, cellIndex = 0;
                        opp_move_status m;

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
    // std::string calcIdx = "";
                        for (const auto& point : vertices) {

                            cellIndex = 0;
                            m = OPP_NEED_MOVE;
    
                            multihop_mover(point, cellIndex, m, 0,0,0);

                            if (m == OPP_MOVE_DONE) { // calcIdx += std::to_string(cellIndex) + " ";
                                mostSuitableCellIndex = std::min(mostSuitableCellIndex, cellIndex);
                            }
                        }    

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "removedCoordinates Writing index %zu | cellIndex %d gbl %d | %s | %s", index, mostSuitableCellIndex, 
    //     (mostSuitableCellIndex < cellSetSizeIncHalo) ? ((int*)global_cell_id_dat->data)[mostSuitableCellIndex]: -1, (mostSuitableCellIndex < cellVolume_dat->set->size) ? "YES" : "NO", calcIdx.c_str());  

                            // TODO : Need an Exclusive MPI_Win_lock for both win_structMeshToCellMapping, win_structMeshToRankMapping windows here

                            // Allow neighbours to write on-behalf of the current rank, to reduce issues
                            if (mostSuitableCellIndex < cellSetSizeIncHalo && 
                                this->structMeshToCellMapping[index] > ((int*)global_cell_id_dat->data)[mostSuitableCellIndex]) {

                                this->structMeshToCellMapping[index] = ((int*)global_cell_id_dat->data)[mostSuitableCellIndex];
#ifdef ENABLE_MPI
                                if (mostSuitableCellIndex < set->size)
                                    this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank
                                else {
                                    
                                    std::map<int, opp_particle_comm_data>* rankVsCommData = nullptr;
                                    for (auto& x : opp_part_comm_neighbour_data) {

                                        if (x.first->cells_set == set) {
                                            rankVsCommData = &(x.second);
                                            break;
                                        }
                                    }

                                    auto it = rankVsCommData->find(mostSuitableCellIndex);
                                    if (it == rankVsCommData->end())
                                        opp_abort(std::string("Cant find comm data of local cell"));

                                    // get the correct neighbour rank
                                    this->structMeshToRankMapping[index] = it->second.cell_residing_rank; 
                                }           
#endif
                            }
                    }
                }
            }

            // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            // structMeshToCellMapping was initialized to MAX_CELL_INDEX, hence reducing to MIN to get the correct cell index
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToRankMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
#endif

            printStructuredMesh("Final");
        }


        //*******************************************************************************
        inline void generateStructMeshToGlobalCellMappings_allSearch(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map) { 

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif            
            if (OP_DEBUG)
                opp_printf("CellMapper", "generateStructMeshToGlobalCellMappings_allSearch global grid dimensions %d %d %d",
                    globalGridDims.x, globalGridDims.y, globalGridDims.z);

            opp_set set = cellVolume_dat->set;

            createStructMeshMappingArrays();

            double lc[N_PER_C];
            int cellSetSizeIncHalo = set->size + set->exec_size + set->nonexec_size;

            // const opp_point& minCoordinate = boundingBox->getLocalMin(); // required for GET_VERT
            const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT

            double x = 0.0, y = 0.0, z = 0.0;
            std::map<size_t, opp_point> removedCoordinates;

            auto all_cell_checker = [&](const opp_point& point, int& cellIndex) {
                
                for (int ci = 0; ci < set->size; ci++) {
                    
                    bool isInside = isPointInCellKernel( 
                                        (const double*)&point, 
                                        lc,
                                        &((double*)cellVolume_dat->data)[ci * cellVolume_dat->dim], 
                                        &((double*)cellDet_dat->data)[ci * cellDet_dat->dim]);
                    
                    if (isInside) {       
                        cellIndex = ci;
                        break;
                    }
                }
            };

            // auto multihop_mover = [&](const opp_point& point, int& cellIndex, opp_move_status& m, int dx, int dy, int dz) {

            //     cellIndex = set->size / 2;

            //     do {
            //         m = getCellIndexKernel(
            //             (const double*)&point, 
            //             &cellIndex,
            //             lc,
            //             &((double*)cellVolume_dat->data)[cellIndex], 
            //             &((double*)cellDet_dat->data)[cellIndex * cellDet_dat->dim], 
            //             &((int*)cellConnectivity_map->map)[cellIndex * cellConnectivity_map->dim]);

            //     } while (m == OPP_NEED_MOVE && cellIndex < cellSetSizeIncHalo);  
            // };

            const opp_point& minGlbCoordinate = boundingBox->getGlobalMin();

            if (cellSetSizeIncHalo > 0) {

                // Step 1 : Get the centroids of the structured mesh cells and try to relate them to unstructured mesh cell indices
                for (int dz = this->localGridStart.z; dz < this->localGridEnd.z; dz++) {
                    
                    z = minGlbCoordinate.z + dz * this->gridSpacing;
                    
                    for (int dy = this->localGridStart.y; dy < this->localGridEnd.y; dy++) {
                        
                        y = minGlbCoordinate.y + dy * this->gridSpacing;
                        
                        for (int dx = this->localGridStart.x; dx < this->localGridEnd.x; dx++) {
                            
                            x = minGlbCoordinate.x + dx * this->gridSpacing;
                            
                            size_t index = (size_t)(dx + dy * globalGridDims.x + dz * globalGridDims.x * globalGridDims.y);
                            
                            const opp_point centroid = getCentroidOfBox(opp_point(x, y ,z));

                            int cellIndex = MAX_CELL_INDEX;

                            all_cell_checker(centroid, cellIndex); // Find in which cell this centroid lies

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "Iterating index %zu | ci %d | setsize %d | %d %d %d | %2.16lE %2.16lE %2.16lE \n| centroid %2.16lE %2.16lE %2.16lE | gbl ci %d", 
    //                                         index, cellIndex, set->size, dx,dy,dz, x,y,z, centroid.x, centroid.y, centroid.z,
    //                                         (cellIndex < set->size) ? ((int*)global_cell_id_dat->data)[cellIndex] : -1);  

                            if (cellIndex == MAX_CELL_INDEX) {
                                removedCoordinates.insert(std::make_pair(index, opp_point(x, y ,z)));
                                continue;
                            }

                            if (cellIndex < set->size) { // write only if the structured cell belong to the current MPI rank

                                this->structMeshToCellMapping[index] = ((int*)global_cell_id_dat->data)[cellIndex];
#ifdef ENABLE_MPI
                                this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank       
#endif
                            }
                        } 
                    }
                }
            }

            // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            // structMeshToCellMapping was initialized to MAX_CELL_INDEX, hence reducing to MIN to get the correct cell index
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToRankMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter

            // The marked structured cells from this rank might be filled by another rank, so if already filled, no need to recalculate from current rank
            for (auto it = removedCoordinates.begin(); it != removedCoordinates.end(); ) {

                size_t removedIndex = it->first;
                
                if (this->structMeshToRankMapping[removedIndex] != MAX_CELL_INDEX) {
                    it = removedCoordinates.erase(it); // This structured index is already written by another rank
                    // opp_printf("CellMapper", "Already written %d to struct index %zu", this->structMeshToRankMapping[removedIndex], removedIndex);
                } else {
                    ++it; 
                }
            }
#endif           

            printStructuredMesh("After centroid mappings");   // Up to this, perfect

            if (cellSetSizeIncHalo > 0) {

                // Step 3 : Iterate over all the NEED_REMOVE points, try to check whether atleast one vertex of the structured mesh can be within 
                //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
                for (auto& p : removedCoordinates) {

                    size_t index = p.first;
                    double& x = p.second.x;
                    double& y = p.second.y;
                    double& z = p.second.z;
                        
                    // still no one has claimed that this cell belongs to it

                    const double gs = gridSpacing;
                    int mostSuitableCellIndex = MAX_CELL_INDEX, mostSuitableGblCellIndex = MAX_CELL_INDEX;

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

// std::string calcIdx = "";
                    for (const auto& point : vertices) {

                        int cellIndex = MAX_CELL_INDEX;

                        all_cell_checker(point, cellIndex);

                        if (cellIndex != MAX_CELL_INDEX && (cellIndex < cellSetSizeIncHalo)) { 

                            int gblCellIndex = ((int*)global_cell_id_dat->data)[cellIndex];
// calcIdx += std::to_string(gblCellIndex) + " ";
                            if (mostSuitableGblCellIndex > gblCellIndex) {
                                mostSuitableGblCellIndex = gblCellIndex;
                                mostSuitableCellIndex = cellIndex;
                            }
                        }
                    }    

    // if (index == 43072) 
    //     opp_printf("ZZZZZ", "removedCoordinates Writing index %zu | cellIndex %d gbl %d | %s | gbl idcs: %s", index, mostSuitableCellIndex, 
    //     mostSuitableGblCellIndex, (mostSuitableCellIndex < cellVolume_dat->set->size) ? "CORE" : "HALO or EX", calcIdx.c_str());  

                    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                    int alreadyAvailGblCellIndex = this->structMeshToCellMapping[index];

                    // Allow neighbours to write on-behalf of the current rank, to reduce issues
                    if (mostSuitableGblCellIndex != MAX_CELL_INDEX && mostSuitableGblCellIndex < alreadyAvailGblCellIndex) {
                        
                        this->structMeshToCellMapping[index] = mostSuitableGblCellIndex;

#ifdef ENABLE_MPI
                        if (mostSuitableCellIndex < set->size) {

                            this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank
                        }
                        else {
                            
                            std::map<int, opp_particle_comm_data>* rankVsCommData = nullptr;
                            for (auto& x : opp_part_comm_neighbour_data) {

                                if (x.first->cells_set == set) {
                                    rankVsCommData = &(x.second);
                                    break;
                                }
                            }

                            auto it = rankVsCommData->find(mostSuitableCellIndex);
                            if (it == rankVsCommData->end())
                                opp_abort(std::string("Cant find comm data of local cell"));

                            // get the correct neighbour rank
                            this->structMeshToRankMapping[index] = it->second.cell_residing_rank; 
                        }           
#endif
                    }

                    CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                    CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
                }
            }

            // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            // structMeshToCellMapping was initialized to MAX_CELL_INDEX, hence reducing to MIN to get the correct cell index
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
            MPI_Allreduce(MPI_IN_PLACE, this->structMeshToRankMapping, globalGridSize, MPI_INT, MPI_MIN, MPI_COMM_WORLD); // comm->comm_inter
#endif

            printStructuredMesh("Final With Global CID");

            for (size_t i = 0; i < globalGridSize; i++) {
                if (this->structMeshToRankMapping[i] == OPP_rank) {
                    const int globalCID = this->structMeshToCellMapping[i];
                    this->structMeshToCellMapping[i] = getLocalCellIndexFromGlobal(globalCID);
                }
            }

            // TODO : For multi node MPI, get this sorted over all ranks

            printStructuredMesh("Final With Local CID");
        }

        //*******************************************************************************
        inline void printStructuredMesh(const std::string msg, bool cellIndices = true) {
            
            if (!OP_DEBUG)
                return;

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (OPP_rank == OPP_ROOT) {

                opp_printf("structMeshToCellMapping", "%s - %s size=%zu", 
                    msg.c_str(), (cellIndices ? "CellIndices" : "MPIRanks"), globalGridSize);

                for (size_t i = 0; i < globalGridSize; i++) {
                    printf("%zu|%d ", i, (cellIndices ? this->structMeshToCellMapping[i] : this->structMeshToRankMapping[i]));
                    if (i % 400 == 0) printf("\n");
                }   
                printf("\n");
            }

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }

        //*******************************************************************************
        inline void createStructMeshMappingArrays() {

            globalGridSize = (size_t)(globalGridDims.x * globalGridDims.y * globalGridDims.z);

#ifdef ENABLE_MPI

            // One per shared memory (node)
            
            // create CELL INDEX array
            {
                const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

                CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                            (void *)&this->structMeshToCellMapping, &this->win_structMeshToCellMapping))

                MPI_Aint allocatedSize = 0;
                int disp = 0;
                CHECK(MPI_Win_shared_query(this->win_structMeshToCellMapping, 0, &allocatedSize, &disp,
                                            (void *)&this->structMeshToCellMapping))

                if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
                    opp_abort(std::string("Pointer to incorrect size in MPI structMeshToCellMapping"));
                }
                if (disp != sizeof(int)) {
                    opp_abort(std::string("Invalid displacement unit in MPI structMeshToCellMapping"));
                }

                if (comm->rank_intra == 0) {
                    for (size_t i = 0; i < globalGridSize; i++)
                        this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
                }

                MPI_Win_fence(0, this->win_structMeshToCellMapping);
            }

            // create RANK array
            {
                const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

                CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                            (void *)&this->structMeshToRankMapping, &this->win_structMeshToRankMapping))

                MPI_Aint allocatedSize = 0;
                int disp = 0;
                CHECK(MPI_Win_shared_query(this->win_structMeshToRankMapping, 0, &allocatedSize, &disp,
                                            (void *)&this->structMeshToRankMapping))

                if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
                    opp_abort(std::string("Pointer to incorrect size in MPI structMeshToRankMapping"));
                }
                if (disp != sizeof(int)) {
                    opp_abort(std::string("Invalid displacement unit in MPI structMeshToRankMapping"));
                }

                if (comm->rank_intra == 0) {
                    for (size_t i = 0; i < globalGridSize; i++)
                        this->structMeshToRankMapping[i] = MAX_CELL_INDEX;
                }

                MPI_Win_fence(0, this->win_structMeshToRankMapping);
            }
#else
            this->structMeshToCellMapping = new int[globalGridSize];
            for (size_t i = 0; i < globalGridSize; i++)
                this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
#endif
        }

        // //*******************************************************************************
        // // This might have double precision issues
        // inline size_t getCellIndexMappingIndex(const opp_point& position) {

        //     int xIndex = static_cast<int>((position.x - this->minGlbCoordinate.x) * this->oneOverGridSpacing);
        //     int yIndex = static_cast<int>((position.y - this->minGlbCoordinate.y) * this->oneOverGridSpacing);
        //     int zIndex = static_cast<int>((position.z - this->minGlbCoordinate.z) * this->oneOverGridSpacing);

        //     return ((size_t)xIndex + (size_t)yIndex * globalGridDims.x + (size_t)zIndex * globalGridDims.x * globalGridDims.y);
        // }

        //*******************************************************************************
        inline size_t getCellIndexMappingIndex(const opp_point& position) {
            // Perform the calculations in higher precision (double)
            double xDiff = position.x - this->minGlbCoordinate.x;
            double yDiff = position.y - this->minGlbCoordinate.y;
            double zDiff = position.z - this->minGlbCoordinate.z;

            xDiff = xDiff * this->oneOverGridSpacing;
            yDiff = yDiff * this->oneOverGridSpacing;
            zDiff = zDiff * this->oneOverGridSpacing;

            // Round to the nearest integer to minimize rounding errors
            const int xIndex = static_cast<int>(xDiff);
            const int yIndex = static_cast<int>(yDiff);
            const int zIndex = static_cast<int>(zDiff);

            // Calculate the cell index mapping index
            return ((size_t)(xIndex) + (size_t)(yIndex * globalGridDims.x) + 
                        (size_t)(zIndex * globalGridDims.x * globalGridDims.y));
        }

        //*******************************************************************************
        // This will contain mappings for halo indices too
        inline void generateGlobalToLocalCellIndexMapping(const opp_dat global_cell_id_dat) {
            
            globalToLocalCellIndexMapping.clear();

            const opp_set cells_set = global_cell_id_dat->set;
            int size_inc_halo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;

            for (int i = 0; i < size_inc_halo; i++) {

                int glbIndex = ((int*)global_cell_id_dat->data)[i];
                
                globalToLocalCellIndexMapping.insert(std::make_pair(glbIndex, i));
            }
        }

        //*******************************************************************************
        inline int getLocalCellIndexFromGlobal(const int globalIndex) {

            if (globalIndex == MAX_CELL_INDEX)
                return MAX_CELL_INDEX;

            auto it = globalToLocalCellIndexMapping.find(globalIndex);
            if (it != globalToLocalCellIndexMapping.end())
                return it->second;
            
            opp_printf("CellMapper", "Error... local cell index not found for global index %d", globalIndex);
            return MAX_INT;
        }

    private:
        const std::shared_ptr<BoundingBox> boundingBox = nullptr;
        const double gridSpacing = 0.0;
        const double oneOverGridSpacing = 0.0;
        const opp_point& minGlbCoordinate;  
        const std::shared_ptr<Comm> comm = nullptr;

        opp_ipoint globalGridDims;
        opp_ipoint localGridStart, localGridEnd;

        std::map<int,int> globalToLocalCellIndexMapping;

        size_t globalGridSize = 0;

        int* structMeshToCellMapping = nullptr;         // This should contain mapping to global cell indices
        int* structMeshToRankMapping = nullptr;
        
        MPI_Win win_structMeshToCellMapping;
        MPI_Win win_structMeshToRankMapping;

        const double ONE_OVER_SIX = (1.0 / 6.0);
        const int N_PER_C = 4;
        const int DET_FIELDS = 4;
        const int NEIGHB_C = 4;
    };
};