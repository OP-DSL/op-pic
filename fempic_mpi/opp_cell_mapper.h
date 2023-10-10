#pragma once

#include "opp_bounding_box.h"
#include "opp_particle_mover_kernel.h"
#include <fstream>

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

            int ax = 0, ay = 0, az = 0;

            { // removed this and added below due to decimal point issues
                // this->globalGridDims.x = std::ceil((maxGblCoordinate.x - minGblCoordinate.x) * oneOverGridSpacing);
                // this->globalGridDims.y = std::ceil((maxGblCoordinate.y - minGblCoordinate.y) * oneOverGridSpacing);
                // this->globalGridDims.z = std::ceil((maxGblCoordinate.z - minGblCoordinate.z) * oneOverGridSpacing);   
            }
            {
                for (double z = minGblCoordinate.z; z < maxGblCoordinate.z; z += this->gridSpacing) az++;
                for (double y = minGblCoordinate.y; y < maxGblCoordinate.y; y += this->gridSpacing) ay++;
                for (double x = minGblCoordinate.x; x < maxGblCoordinate.x; x += this->gridSpacing) ax++; 
                this->globalGridDims.x = ax + 1;
                this->globalGridDims.y = ay + 1;
                this->globalGridDims.z = az + 1; 
            }

            if (OPP_rank == OPP_ROOT)
                opp_printf("CellMapper", "Global Grid Size - [%d %d %d] gridSpacing [%2.10lE]", 
                    this->globalGridDims.x, this->globalGridDims.y, this->globalGridDims.z, this->gridSpacing); 
            
            const opp_point& minLocalCoordinate = boundingBox->getLocalMin();
            const opp_point& maxLocalCoordinate = boundingBox->getLocalMax();

            // Find the local ranks grid start indices
            ax = 0; ay = 0; az = 0;
            for (double z = minGblCoordinate.z; (z < minLocalCoordinate.z); z += this->gridSpacing) az++;
            for (double y = minGblCoordinate.y; (y < minLocalCoordinate.y); y += this->gridSpacing) ay++; 
            for (double x = minGblCoordinate.x; (x < minLocalCoordinate.x); x += this->gridSpacing) ax++; 
            this->localGridStart.x = ax == 0 ? ax : (ax - 1);
            this->localGridStart.y = ay == 0 ? ay : (ay - 1);
            this->localGridStart.z = az == 0 ? az : (az - 1);         

            // Find the local ranks grid end indices
            ax = 0; ay = 0; az = 0;
            for (double z = minGblCoordinate.z; (z <= maxLocalCoordinate.z); z += this->gridSpacing) az++; 
            for (double y = minGblCoordinate.y; (y <= maxLocalCoordinate.y); y += this->gridSpacing) ay++; 
            for (double x = minGblCoordinate.x; (x <= maxLocalCoordinate.x); x += this->gridSpacing) ax++; 
            this->localGridEnd.x = this->globalGridDims.x == ax ? ax : (ax + 1);
            this->localGridEnd.y = this->globalGridDims.y == ay ? ay : (ay + 1);
            this->localGridEnd.z = this->globalGridDims.z == az ? az : (az + 1);     

            if (OP_DEBUG)
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
        inline int findClosestLocalCellIndex(const size_t& structCellIdx) { 
            
            // if (OP_DEBUG) 
            {
                if (structCellIdx >= globalGridSize) {
                    opp_printf("findClosestLocalCellIndex", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                        structCellIdx, globalGridSize);
                    return MAX_CELL_INDEX;
                }
            }
            
            return this->structMeshToCellMapping[structCellIdx];
        }

        //*******************************************************************************
        // Returns the rank of the cell
        inline int findClosestCellRank(const size_t& structCellIdx) { 

#ifdef ENABLE_MPI 
            // if (OP_DEBUG) 
            {
                if (structCellIdx >= globalGridSize) {
                    opp_printf("findClosestCellRank", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                        structCellIdx, globalGridSize);
                    return MAX_CELL_INDEX;
                }
            }

            return this->structMeshToRankMapping[structCellIdx];
#else
            return OPP_rank;
#endif
        }

// Be careful when implementing MPI
// 1: structMeshToCellMapping should contain global indices - DONE
// 2: it should be populated one copy per shared memory (need to handle overlappings)
// 3: once populated per shared memory, need to share the details with inter nodes (can use inter node comm root) - DONE, CHECK

        //*******************************************************************************
        inline void generateStructMeshToGlobalCellMappings(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map) { 

            generateStructMeshToGlobalCellMappings_allSearch(cellVolume_dat, cellDet_dat, global_cell_id_dat, cellConnectivity_map);
        }

        //*******************************************************************************
        inline void reduceInterNodeMappings(int callID) {

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);

            if (comm->rank_intra == 0) { // comm->size_intra > 1 && 

                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                if (OP_DEBUG) 
                    opp_printf("CellMapper", "reduceInterNodeMappings Start size_inter %d", comm->size_inter);

                std::string sendLog = "Send Counts: ", recvLog = "Recv Counts: ";

                int totalSendCount = 0, totalRecvCount = 0;
                std::vector<int> sendCounts(comm->size_inter);
                std::vector<int> sendDisplacements(comm->size_inter);
                std::vector<int> recvCounts(comm->size_inter);
                std::vector<int> recvDisplacements(comm->size_inter);

                for (int i = 0; i < comm->size_inter; ++i) {
                    const int remainder = (int)(this->globalGridSize % comm->size_inter);
                    sendCounts[i] = (this->globalGridSize / comm->size_inter) + (i < remainder ? 1 : 0);
                    sendDisplacements[i] = totalSendCount;
                    totalSendCount += sendCounts[i];

                    sendLog += std::to_string(sendCounts[i]) + "|Disp:" + std::to_string(sendDisplacements[i]) + " ";
                }

                for (int i = 0; i < comm->size_inter; ++i) {
                    recvCounts[i] = sendCounts[comm->rank_inter];
                    recvDisplacements[i] = totalRecvCount;
                    totalRecvCount += recvCounts[i];

                    recvLog += std::to_string(recvCounts[i]) + "|Disp:" + std::to_string(recvDisplacements[i]) + " ";
                }

                // if (OP_DEBUG) 
                {
                    opp_printf("CellMapper", "reduceInterNodeMappings totalSendCount %d : %s", totalSendCount, sendLog.c_str());
                    opp_printf("CellMapper", "reduceInterNodeMappings totalRecvCount %d : %s", totalRecvCount, recvLog.c_str());
                }

                std::vector<int> cellMappingsRecv(totalRecvCount, MAX_CELL_INDEX-1);
                std::vector<int> ranksRecv(totalRecvCount, MAX_CELL_INDEX-1);

                const int INT_COUNT_PER_MSG = 10000;

                std::vector<MPI_Request> sendRequests;
                for (int rank = 0; rank < comm->size_inter; rank++) 
                {
                    if (rank == comm->rank_inter) continue;

                    int totalSendCount = sendCounts[rank];
                    int sendCount = 0, alreadySentCount = 0;

                    int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

                    for (int i = 0; i < blocks; i++) {
                        
                        if (i != blocks-1) 
                            sendCount = INT_COUNT_PER_MSG;
                        else 
                            sendCount = totalSendCount - alreadySentCount;

                        // opp_printf("SEND", "blocks %d|%d disp %d sendCounts %d to rank %d | %d %d %d %d %d %d %d %d %d %d", 
                        //     i, blocks,
                        //     sendDisplacements[rank] + i * INT_COUNT_PER_MSG, sendCount, inter_ranks[rank],
                        //     sb[0], sb[1], sb[2], sb[3], sb[4], sb[5], sb[6], sb[7], sb[8], sb[9]);

                        sendRequests.emplace_back(MPI_Request());
                        MPI_Isend(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                                rank, (10000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                        sendRequests.emplace_back(MPI_Request());
                        MPI_Isend(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                                rank, (20000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                        alreadySentCount += sendCount;
                    }
                }

                std::vector<MPI_Request> recvRequests;
                for (int rank = 0; rank < comm->size_inter; rank++) 
                {
                    if (rank == comm->rank_inter) continue;

                    int totalRecvCount = recvCounts[rank];
                    int recvCount = 0, alreadyRecvCount = 0;

                    int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

                    for (int i = 0; i < blocks; i++) {

                        if (i != blocks-1) 
                            recvCount = INT_COUNT_PER_MSG;
                        else 
                            recvCount = totalRecvCount - alreadyRecvCount;

                        // opp_printf("RECV", "blocks %d|%d disp %d recvCounts %d from rank %d", i, blocks,
                        //     recvDisplacements[rank] + i * INT_COUNT_PER_MSG, recvCount, inter_ranks[rank]);

                        recvRequests.emplace_back(MPI_Request());
                        MPI_Irecv(&(ranksRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                                rank, (10000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                        recvRequests.emplace_back(MPI_Request());
                        MPI_Irecv(&(cellMappingsRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                                rank, (20000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                        alreadyRecvCount += recvCount;
                    }
                }

                // opp_printf("DISPLACEMENTS", "send %d recv %d", sendDisplacements[comm->rank_inter], recvDisplacements[comm->rank_inter]);

                // Copy own data
                for (int i = 0; i < recvCounts[comm->rank_inter]; i++) {

                    ranksRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i];
                    cellMappingsRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i];
                }

                // opp_printf("CellMapper", "Going for MPI_Wait");
                MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
                MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
                MPI_Barrier(comm->comm_inter);
                // opp_printf("CellMapper", "We are proceeeeeeding");

                // MPI_Alltoallv(this->structMeshToCellMapping, sendCounts.data(), sendDisplacements.data(), MPI_INT, 
                //     cellMappingsRecv.data(), recvCounts.data(), recvDisplacements.data(), MPI_INT, comm->comm_inter);
                // MPI_Alltoallv(this->structMeshToRankMapping, sendCounts.data(), sendDisplacements.data(), MPI_INT, 
                //     ranksRecv.data(), recvCounts.data(), recvDisplacements.data(), MPI_INT, comm->comm_inter);
                // MPI_Barrier(comm->comm_inter);

                // printStructuredMesh(std::string("RECV_MAPPING") + std::to_string(callID), cellMappingsRecv.data(), totalRecvCount);
                // printStructuredMesh(std::string("RECV_RANKS") + std::to_string(callID), ranksRecv.data(), totalRecvCount);
                // MPI_Barrier(comm->comm_inter);

                const size_t recvCount = recvCounts[0];
                for (size_t i = 0; i < recvCount; i++) {

                    int cellIndex = MAX_CELL_INDEX;
                    for (int r = 0; r < comm->size_inter; r++) {

                        if (cellIndex > cellMappingsRecv[i + r * recvCount]) {

                            cellIndex = cellMappingsRecv[i + r * recvCount];
                            cellMappingsRecv[i] = cellIndex;
                            ranksRecv[i] = ranksRecv[i + r * recvCount];
                        }
                    }
                }

                MPI_Barrier(comm->comm_inter);

                // printStructuredMesh(std::string("MAPPING_COMP") + std::to_string(callID), ranksRecv.data(), recvCount);
                // MPI_Barrier(comm->comm_inter);

                sendRequests.clear();
                for (int rank = 0; rank < comm->size_inter; rank++) 
                {
                    if (rank == comm->rank_inter) continue;

                    int totalSendCount = recvCount;
                    int sendCount = 0, alreadySentCount = 0;

                    int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

                    for (int i = 0; i < blocks; i++) {
                        
                        if (i != blocks-1) 
                            sendCount = INT_COUNT_PER_MSG;
                        else 
                            sendCount = totalSendCount - alreadySentCount;

                        sendRequests.emplace_back(MPI_Request());
                        MPI_Isend(&(cellMappingsRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (30000 + i), 
                                comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                        sendRequests.emplace_back(MPI_Request());
                        MPI_Isend(&(ranksRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (40000 + i), 
                                comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                        alreadySentCount += sendCount;
                    }
                }

                // Since we are about to write to data mapped by windows, release the lock here
                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));

                recvRequests.clear();
                for (int rank = 0; rank < comm->size_inter; rank++) 
                {
                    if (rank == comm->rank_inter) continue;

                    int totalRecvCount = sendCounts[rank];
                    int recvCount2 = 0, alreadyRecvCount = 0;

                    int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

                    for (int i = 0; i < blocks; i++) {

                        if (i != blocks-1) 
                            recvCount2 = INT_COUNT_PER_MSG;
                        else 
                            recvCount2 = totalRecvCount - alreadyRecvCount;

                        recvRequests.emplace_back(MPI_Request());
                        MPI_Irecv(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                                rank, (30000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                        recvRequests.emplace_back(MPI_Request());
                        MPI_Irecv(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                                rank, (40000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                        alreadyRecvCount += recvCount2;
                    }
                }

                // Copy own data
                for (int i = 0; i < sendCounts[comm->rank_inter]; i++) {

                    this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i] = ranksRecv[i];
                    this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i] = cellMappingsRecv[i];
                }

                // opp_printf("CellMapper", "Going for MPI_Wait 2");
                MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
                MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
                MPI_Barrier(comm->comm_inter);
                // opp_printf("CellMapper", "We are proceeeeeeding 2");

                // MPI_Allgatherv(cellMappingsRecv.data(), recvCount, MPI_INT, this->structMeshToCellMapping, 
                //     sendCounts.data(), sendDisplacements.data(), MPI_INT, comm->comm_inter);         
                // MPI_Allgatherv(ranksRecv.data(), recvCount, MPI_INT, this->structMeshToRankMapping, 
                //     sendCounts.data(), sendDisplacements.data(), MPI_INT, comm->comm_inter);
                // CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                // CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
                
                if (OP_DEBUG) 
                    opp_printf("CellMapper", "reduceInterNodeMappings END");
            }

            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
        
        //*******************************************************************************
        inline void convertToLocalMappings() {

            if (OP_DEBUG) 
            // if (OPP_rank == 0)
                opp_printf("CellMapper", "convertToLocalMappings Start");

            if (!globalCellMappingInitialized) {
                opp_abort("Global Cell Mappings not initialized (@convertToLocalMappings)");
            }

#ifdef ENABLE_MPI
            if (comm->rank_intra == 0) {
                
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));

                for (size_t i = 0; i < globalGridSize; i++) {
                    
                    if (this->structMeshToCellMapping[i] != MAX_CELL_INDEX)
                        this->structMeshToCellMapping[i] = (-1 * this->structMeshToCellMapping[i]);
                }

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
            }

            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
#endif

            // bool errorFound = false;

            for (size_t i = 0; i < globalGridSize; i++) {
                
                if (this->structMeshToRankMapping[i] == OPP_rank) {
                    
                    const int globalCID = (-1 * this->structMeshToCellMapping[i]);
                    
                    if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {
                        
                        const int localCID = getLocalCellIndexFromGlobal(globalCID);
                        
                        if (localCID != MAX_CELL_INDEX) {
                            this->structMeshToCellMapping[i] = localCID;
                        }
                        else { // if (OP_DEBUG) {
                            // errorFound = true;
                            opp_printf("CellMapper", "Error at convertToLocalMappings : structMeshToCellMapping at index %d is invalid [gcid:%d]", 
                                i, this->structMeshToCellMapping[i]);
                        }
                    }
                }
            }

            // if (errorFound) {
            //     opp_abort("Negative structMeshToCellMappings found (@convertToLocalMappings)");
            // }

#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);

            if (comm->rank_intra == 0) {

                MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MAX, comm->comm_inter);
            }

            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (OP_DEBUG) 
            // if (OPP_rank == 0)
                opp_printf("CellMapper", "convertToLocalMappings END");
        }

        //*******************************************************************************
        inline void generateStructMeshToGlobalCellMappings_allSearch(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map) { 

            // if (OP_DEBUG) 
            if (OPP_rank == 0)            
                opp_printf("CellMapper", "generateStructMeshToGlobalCellMappings_allSearch start");

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

                            if (cellIndex == MAX_CELL_INDEX) {
                                removedCoordinates.insert(std::make_pair(index, opp_point(x, y ,z)));
                                continue;
                            }

                            if (cellIndex < set->size) { // write only if the structured cell belong to the current MPI rank

                                MPI_Put(&((int*)global_cell_id_dat->data)[cellIndex], 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToCellMapping);
                                // this->structMeshToCellMapping[index] = ((int*)global_cell_id_dat->data)[cellIndex];
#ifdef ENABLE_MPI
                                MPI_Put(&(comm->rank_parent), 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToRankMapping);
                                // this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank       
#endif
                            }
                        } 
                    }
                }
            }

            // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            if (comm->rank_intra == 0) {

                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));
                
                printStructuredMesh(std::string("Core_Bef_M"), structMeshToCellMapping, globalGridSize);
                printStructuredMesh(std::string("Core_Bef_R"), structMeshToRankMapping, globalGridSize);

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
            }

            MPI_Barrier(MPI_COMM_WORLD);

            reduceInterNodeMappings(1);

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            MPI_Barrier(MPI_COMM_WORLD);

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

            MPI_Barrier(MPI_COMM_WORLD);

            if (comm->rank_intra == 0) {
                
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                printStructuredMesh(std::string("Core_Aft_M"), structMeshToCellMapping, globalGridSize);
                printStructuredMesh(std::string("Core_Aft_R"), structMeshToRankMapping, globalGridSize);

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 
            
#endif           

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

                    for (const auto& point : vertices) {

                        int cellIndex = MAX_CELL_INDEX;

                        all_cell_checker(point, cellIndex);

                        if (cellIndex != MAX_CELL_INDEX && (cellIndex < cellSetSizeIncHalo)) { 

                            int gblCellIndex = ((int*)global_cell_id_dat->data)[cellIndex];

                            if (mostSuitableGblCellIndex > gblCellIndex) {
                                mostSuitableGblCellIndex = gblCellIndex;
                                mostSuitableCellIndex = cellIndex;
                            }
                        }
                    }    

                    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                    int alreadyAvailGblCellIndex = this->structMeshToCellMapping[index];

                    // Allow neighbours to write on-behalf of the current rank, to reduce issues
                    if (mostSuitableGblCellIndex != MAX_CELL_INDEX && mostSuitableGblCellIndex < alreadyAvailGblCellIndex && 
                        (mostSuitableCellIndex < set->size)) {

                        this->structMeshToCellMapping[index] = mostSuitableGblCellIndex;
#ifdef ENABLE_MPI
                        this->structMeshToRankMapping[index] = comm->rank_parent; // Global MPI rank          
#endif
                    }

                    CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                    CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
                }
            }

            // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
#ifdef ENABLE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            if (comm->rank_intra == 0) {
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                printStructuredMesh(std::string("Edge_Bef_M"), structMeshToCellMapping, globalGridSize);
                printStructuredMesh(std::string("Edge_Bef_R"), structMeshToRankMapping, globalGridSize);

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
            }

            MPI_Barrier(MPI_COMM_WORLD);

            reduceInterNodeMappings(2);

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            if (comm->rank_intra == 0) {
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                printStructuredMesh(std::string("Edge_Aft_M"), structMeshToCellMapping, globalGridSize);
                printStructuredMesh(std::string("Edge_Aft_R"), structMeshToRankMapping, globalGridSize);

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            convertToLocalMappings();

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Win_fence(0, this->win_structMeshToCellMapping); 
            MPI_Win_fence(0, this->win_structMeshToRankMapping); 

            if (comm->rank_intra == 0) {
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

                printStructuredMesh(std::string("Aft_ConvLocal_M"), structMeshToCellMapping, globalGridSize);
                printStructuredMesh(std::string("Aft_ConvLocal_R"), structMeshToRankMapping, globalGridSize);

                CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
                CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
            }
#endif

            // if (OP_DEBUG) 
            if (OPP_rank == 0)
                opp_printf("CellMapper", "generateStructMeshToGlobalCellMappings_allSearch end");
        }

        //*******************************************************************************
        inline void printStructuredMesh(const std::string msg, int *array, size_t size, bool printToFile = true) {
            
            // if (!OP_DEBUG)
            //     return;

            if (!printToFile) 
            {
                opp_printf("structMeshToCellMapping", "%s - size=%zu", msg.c_str(), size);

                for (size_t i = 0; i < size; i++) {
                    printf("%zu|%d ", i, array[i]);
                    if (i % 50 == 0) printf("\n");
                }   
                printf("\n");
            }
            else
            {
                const std::string file_name = std::string("files/struct_com") + std::to_string(OPP_comm_size) + "_r" + 
                                                std::to_string(OPP_rank) + "_" + msg; 

                std::ofstream outFile(file_name);

                if (!outFile.is_open()) {
                    opp_printf("printStructuredMesh", "can't open file %s\n", file_name.c_str());
                    opp_abort("printStructuredMesh - can't open file");
                }

                for (size_t i = 0; i < size; i++) {

                    outFile << i << "|";

                    int value = array[i];
                    if (value != MAX_CELL_INDEX) {
                        outFile << value << " ";
                    }
                    else {
                        outFile << "X ";
                    }

                    if ((i+1) % 50 == 0) 
                        outFile << "\n";
                }   
                outFile << "\n";
                outFile.close();
            }
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

        //*******************************************************************************
        inline size_t getCellIndexMappingIndex(const opp_point& position) {

            // int xIndex = static_cast<int>((position.x - this->minGlbCoordinate.x) * this->oneOverGridSpacing);
            // int yIndex = static_cast<int>((position.y - this->minGlbCoordinate.y) * this->oneOverGridSpacing);
            // int zIndex = static_cast<int>((position.z - this->minGlbCoordinate.z) * this->oneOverGridSpacing);

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
            size_t index = ((size_t)(xIndex) + (size_t)(yIndex * globalGridDims.x) + 
                        (size_t)(zIndex * globalGridDims.x * globalGridDims.y));

            if (index >= globalGridSize) {
                opp_printf("CellMapper", "Ooi ooi index %zu globalGridSize %zu", index, globalGridSize);
            }

            return index;
        }

        //*******************************************************************************
        // This will contain mappings for halo indices too
        inline void generateGlobalToLocalCellIndexMapping(const opp_dat global_cell_id_dat) {

            // if (OP_DEBUG) 
            if (OPP_rank == 0)            
                opp_printf("CellMapper", "generateGlobalToLocalCellIndexMapping start");
                
            globalToLocalCellIndexMapping.clear();

            const opp_set cells_set = global_cell_id_dat->set;
            int size_inc_halo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;

            for (int i = 0; i < size_inc_halo; i++) {

                int glbIndex = ((int*)global_cell_id_dat->data)[i];
                
                globalToLocalCellIndexMapping.insert(std::make_pair(glbIndex, i));
            }
            
            globalCellMappingInitialized = true;

            // if (OP_DEBUG) 
            if (OPP_rank == 0)
                opp_printf("CellMapper", "generateGlobalToLocalCellIndexMapping end");
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

        //*******************************************************************************
        inline void destroyGlobalToLocalCellIndexMapping() {
            
            globalToLocalCellIndexMapping.clear();

            if (OPP_rank == 0)
                opp_printf("CellMapper", "destroyGlobalToLocalCellIndexMapping");
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
        bool globalCellMappingInitialized = false;

        size_t globalGridSize = 0;
        
        int* structMeshToCellMapping = nullptr;         // This contain mapping to local cell indices
        int* structMeshToRankMapping = nullptr;
        
        MPI_Win win_structMeshToCellMapping;
        MPI_Win win_structMeshToRankMapping;

        const double ONE_OVER_SIX = (1.0 / 6.0);
        const int N_PER_C = 4;
        const int DET_FIELDS = 4;
        const int NEIGHB_C = 4;
    };
};