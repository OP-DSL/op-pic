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

#include <opp_lib_core.h>

#define BOUNDING_TOLERENCE 1e-12

//*******************************************************************************
using namespace opp;

std::shared_ptr<BoundingBox> boundingBox;
std::shared_ptr<CellMapper> cellMapper;
std::shared_ptr<Comm> comm;
std::unique_ptr<GlobalParticleMover> globalMover;
bool useGlobalMove = false;

//*******************************************************************************
BoundingBox::BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate) :
    dim(dim) {

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(0);

    if (OPP_DBG)
        opp_printf("Local BoundingBox [provided]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);

}

// For now, only 3D is implemented
//*******************************************************************************
BoundingBox::BoundingBox(const opp_dat node_pos_dat, int dim) :
    dim(dim) {
    
    if (dim != 3 && dim != 2) {
        opp_abort(std::string("For now, only 2D/3D BoundingBox is implemented"));
    }

    const double* node_pos_data = (const double*)node_pos_dat->data;
    const int node_count = node_pos_dat->set->size + node_pos_dat->set->exec_size + 
                                node_pos_dat->set->nonexec_size;;

{
    opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
    opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

    if (dim == 2) {
        minCoordinate.z = 0;
        maxCoordinate.z = 0;
    }

    // make the bounding box even over the halo regions
    for (int i = 0; i < node_pos_dat->set->size; i++) {
        minCoordinate.x = std::min(node_pos_data[i * dim + 0], minCoordinate.x);
        maxCoordinate.x = std::max(node_pos_data[i * dim + 0], maxCoordinate.x);
        minCoordinate.y = std::min(node_pos_data[i * dim + 1], minCoordinate.y);
        maxCoordinate.y = std::max(node_pos_data[i * dim + 1], maxCoordinate.y);
        if (dim == 3) {
            minCoordinate.z = std::min(node_pos_data[i * dim + 2], minCoordinate.z);
            maxCoordinate.z = std::max(node_pos_data[i * dim + 2], maxCoordinate.z);
        }
    }

        opp_printf("Local BoundingBox [computed]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            minCoordinate.x, minCoordinate.y, minCoordinate.z, 
            minCoordinate.x, minCoordinate.y, minCoordinate.z);
}

    opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
    opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

    if (dim == 2) {
        minCoordinate.z = 0;
        maxCoordinate.z = 0;
    }

    // make the bounding box even over the halo regions
    for (int i = 0; i < node_count; i++) {
        minCoordinate.x = std::min(node_pos_data[i * dim + 0], minCoordinate.x);
        maxCoordinate.x = std::max(node_pos_data[i * dim + 0], maxCoordinate.x);
        minCoordinate.y = std::min(node_pos_data[i * dim + 1], minCoordinate.y);
        maxCoordinate.y = std::max(node_pos_data[i * dim + 1], maxCoordinate.y);
        if (dim == 3) {
            minCoordinate.z = std::min(node_pos_data[i * dim + 2], minCoordinate.z);
            maxCoordinate.z = std::max(node_pos_data[i * dim + 2], maxCoordinate.z);
        }
    }

    if (node_count != 0) {
        minCoordinate.x -= BOUNDING_TOLERENCE;
        maxCoordinate.x += BOUNDING_TOLERENCE;
        minCoordinate.y -= BOUNDING_TOLERENCE;
        maxCoordinate.y += BOUNDING_TOLERENCE;
        if (dim == 3) {
            minCoordinate.z -= BOUNDING_TOLERENCE;
            maxCoordinate.z += BOUNDING_TOLERENCE;
        }
    }

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(node_count);

    if (OPP_DBG)
        opp_printf("Local BoundingBox [computed]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);
}

//*******************************************************************************
void BoundingBox::generateGlobalBoundingBox(int count) {

#ifdef USE_MPI

    const double* localMin = reinterpret_cast<const double*>(&(this->boundingBox[0]));
    double* globalMin = reinterpret_cast<double*>(&(this->globalBoundingBox[0]));    
    MPI_Allreduce(localMin, globalMin, dim, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
    double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
    MPI_Allreduce(localMax, globalMax, dim, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // This is to avoid min and max corrdinates to not have MAX_REAL and MIN_REAL when current rank has no work
    if (count == 0) {
        this->boundingBox[0].x = this->globalBoundingBox[0].x; 
        this->boundingBox[0].y = this->globalBoundingBox[0].y; 
        this->boundingBox[0].z = this->globalBoundingBox[0].z; 

        this->boundingBox[1].x = this->globalBoundingBox[0].x;
        this->boundingBox[1].y = this->globalBoundingBox[0].y;
        this->boundingBox[1].z = this->globalBoundingBox[0].z;
    }

    opp_printf("Local BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
        this->boundingBox[0].x, this->boundingBox[0].y, this->boundingBox[0].z, 
        this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);

    if (OPP_rank == OPP_ROOT)
        opp_printf("Global BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->globalBoundingBox[0].x, this->globalBoundingBox[0].y, this->globalBoundingBox[0].z, 
            this->globalBoundingBox[1].x, this->globalBoundingBox[1].y, this->globalBoundingBox[1].z);
#else
    this->globalBoundingBox = this->boundingBox;
#endif
}

//*******************************************************************************
BoundingBox::~BoundingBox() { };

//*******************************************************************************
const opp_point& BoundingBox::getLocalMin() const {

    return this->boundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getLocalMax() const {

    return this->boundingBox[1];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMin() const {

    return this->globalBoundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMax() const {
    
    return this->globalBoundingBox[1];
}

//*******************************************************************************
bool BoundingBox::isCoordinateInBoundingBox(const opp_point& point) { 

    if (this->boundingBox[0].x > point.x || this->boundingBox[1].x < point.x) 
        return false;
    else if (this->boundingBox[0].y > point.y || this->boundingBox[1].y < point.y) 
        return false;
    else if (this->boundingBox[0].z > point.z || this->boundingBox[1].z < point.z) 
        return false;
    
    return true;
}

//*******************************************************************************
bool BoundingBox::isCoordinateInGlobalBoundingBox(const opp_point& point) { 

    if (this->globalBoundingBox[0].x > point.x || this->globalBoundingBox[1].x < point.x) 
        return false;
    else if (this->globalBoundingBox[0].y > point.y || this->globalBoundingBox[1].y < point.y) 
        return false;
    else if (this->globalBoundingBox[0].z > point.z || this->globalBoundingBox[1].z < point.z) 
        return false;
    
    return true;
}


//*******************************************************************************
// This will contain mappings for halo indices too
GlobalToLocalCellIndexMapper::GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat) {

    // if (OPP_DBG) 
    if (OPP_rank == 0)            
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping start");
        
    globalToLocalCellIndexMap.clear();

    const opp_set cells_set = global_cell_id_dat->set;
    int size_inc_halo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;

    for (int i = 0; i < size_inc_halo; i++) {

        int glbIndex = ((int*)global_cell_id_dat->data)[i];
        
        globalToLocalCellIndexMap.insert(std::make_pair(glbIndex, i));
    }
    
    // if (OPP_DBG) 
    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping end");
}

//*******************************************************************************
GlobalToLocalCellIndexMapper::~GlobalToLocalCellIndexMapper() {
    
    globalToLocalCellIndexMap.clear();

    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "destroyGlobalToLocalCellIndexMapping");
}

//*******************************************************************************
int GlobalToLocalCellIndexMapper::map(const int globalIndex) {

    if (globalIndex == MAX_CELL_INDEX)
        return MAX_CELL_INDEX;

    auto it = globalToLocalCellIndexMap.find(globalIndex);
    if (it != globalToLocalCellIndexMap.end())
        return it->second;
    
    opp_printf("GlobalToLocalCellIndexMapper", "Error... local cell index not found for global index %d", globalIndex);
    return MAX_INT;
}


//*******************************************************************************
CellMapper::CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, const std::shared_ptr<Comm> comm) 
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

    if (OPP_DBG)
        opp_printf("CellMapper", "Local Grid - Min[%d %d %d] Max[%d %d %d]", 
            this->localGridStart.x, this->localGridStart.y, this->localGridStart.z, 
            this->localGridEnd.x, this->localGridEnd.y, this->localGridEnd.z); 
}

//*******************************************************************************
CellMapper::~CellMapper() { 

#ifdef USE_MPI
    CHECK(MPI_Win_free(&this->win_structMeshToCellMapping))
    this->structMeshToCellMapping = nullptr;

    CHECK(MPI_Win_free(&this->win_structMeshToRankMapping))
    this->structMeshToRankMapping = nullptr;           
#else
    delete[] this->structMeshToCellMapping;
#endif
};

//*******************************************************************************
opp_point CellMapper::getCentroidOfBox(const opp_point& coordinate) { 

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
size_t CellMapper::findStructuredCellIndex3D(const opp_point& position) { 

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
        // opp_printf("CellMapper", "Error index %zu generated is larger than globalGridSize %zu", 
        //     index, globalGridSize);
        return MAX_CELL_INDEX;
    }

    return index;
}

//*******************************************************************************
// Returns the global cell index
size_t CellMapper::findStructuredCellIndex2D(const opp_point& position) { 

    // Perform the calculations in higher precision (double)
    double xDiff = position.x - this->minGlbCoordinate.x;
    double yDiff = position.y - this->minGlbCoordinate.y;

    xDiff = xDiff * this->oneOverGridSpacing;
    yDiff = yDiff * this->oneOverGridSpacing;

    // Round to the nearest integer to minimize rounding errors
    const int xIndex = static_cast<int>(xDiff);
    const int yIndex = static_cast<int>(yDiff);

    // Calculate the cell index mapping index
    size_t index = ((size_t)(xIndex) + (size_t)(yIndex * globalGridDims.x));

    if (index >= globalGridSize) {
        // opp_printf("CellMapper", "Error index %zu generated is larger than globalGridSize %zu", 
        //     index, globalGridSize);
        return MAX_CELL_INDEX;
    }

    return index;
}

//*******************************************************************************
// Returns the global cell index
int CellMapper::findClosestCellIndex(const size_t& structCellIdx) { 
    
    if (OPP_DBG) 
    {
        if (structCellIdx >= globalGridSize) {
            opp_printf("findClosestCellIndex", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                structCellIdx, globalGridSize);
            return MAX_CELL_INDEX;
        }
    }
    
    return this->structMeshToCellMapping[structCellIdx];
}

//*******************************************************************************
// Returns the rank of the cell
int CellMapper::findClosestCellRank(const size_t& structCellIdx) { 

#ifdef USE_MPI 
    if (OPP_DBG) 
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

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) {

#ifdef USE_MPI
    waitBarrier();

    if (comm->rank_intra == 0) { // comm->size_intra > 1 && 

        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

        if (OPP_DBG) 
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

            if (OPP_DBG) 
                sendLog += std::to_string(sendCounts[i]) + "|Disp:" + std::to_string(sendDisplacements[i]) + " ";
        }

        for (int i = 0; i < comm->size_inter; ++i) {
            recvCounts[i] = sendCounts[comm->rank_inter];
            recvDisplacements[i] = totalRecvCount;
            totalRecvCount += recvCounts[i];

            if (OPP_DBG) 
                recvLog += std::to_string(recvCounts[i]) + "|Disp:" + std::to_string(recvDisplacements[i]) + " ";
        }

        if (OPP_DBG) 
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

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);

        // printStructuredMesh(std::string("RECV_MAPPING") + std::to_string(callID), cellMappingsRecv.data(), totalRecvCount);
        // printStructuredMesh(std::string("RECV_RANKS") + std::to_string(callID), ranksRecv.data(), totalRecvCount);
        // MPI_Barrier(comm->comm_inter);

        const size_t recvCount = recvCounts[0];
        for (size_t i = 0; i < recvCount; i++) { // reduce to get common mapping on all inter node ranks

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

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);
        
        if (OPP_DBG) 
            opp_printf("CellMapper", "reduceInterNodeMappings END");
    }

    waitBarrier();
#endif
}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat);

#ifdef USE_MPI
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

    for (size_t i = 0; i < globalGridSize; i++) {
        
        if (this->structMeshToRankMapping[i] == OPP_rank) {
            
            const int globalCID = (-1 * this->structMeshToCellMapping[i]);
            
            if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {
                
                const int localCID = globalToLocalCellIndexMapper.map(globalCID);
                
                if (localCID != MAX_CELL_INDEX) {
                    this->structMeshToCellMapping[i] = localCID;
                }
                else {
                    opp_printf("CellMapper", "Error at convertToLocalMappings : structMeshToCellMapping at index %d is invalid [gcid:%d]", 
                        i, this->structMeshToCellMapping[i]);
                }
            }
        }
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);

    if (comm->rank_intra == 0) {

        MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MAX, comm->comm_inter);
    }

    waitBarrier();
#endif
    if (OPP_DBG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::enrichStructuredMesh(const int index, const int cell_index, const int rank) {

#ifdef USE_MPI
    MPI_Put(&cell_index, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToCellMapping);
    MPI_Put(&rank, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToRankMapping);
#else
    this->structMeshToCellMapping[index] = cell_index;
#endif
}

//*******************************************************************************
void CellMapper::printStructuredMesh(const std::string msg, int *array, size_t size, bool printToFile) {
    
    // if (!OPP_DBG)
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
void CellMapper::createStructMeshMappingArrays() {

    globalGridSize = (size_t)(globalGridDims.x * globalGridDims.y * globalGridDims.z);

#ifdef USE_MPI

    MPI_Barrier(MPI_COMM_WORLD);

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
void CellMapper::waitBarrier() {

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, this->win_structMeshToCellMapping); 
    MPI_Win_fence(0, this->win_structMeshToRankMapping); 
#endif
}

//*******************************************************************************
void CellMapper::lockWindows() {

#ifdef USE_MPI
    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));
#endif
}


void CellMapper::unlockWindows() {

#ifdef USE_MPI
    CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
#endif
}

