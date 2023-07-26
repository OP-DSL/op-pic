#pragma once

#include "opp_bounding_box.h"
#include "opp_particle_mover_kernel.h"

#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

#define ASSIGN_CENTROID_TO_DIM(K)                                   \
    if (coordinate.K + this->gridSpacing <= maxCoordinate.K) {      \
        centroid.K = coordinate.K + this->gridSpacing / 2;;         \
    }                                                               \
    else {                                                          \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5;        \
    }                                                               \

#ifdef ENABLE_MPI
    #include "opp_comm.h"
#else
    #define Comm void
#endif



namespace opp {

    class CellMapper {
    
    public:

        //*******************************************************************************
        CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, Comm* comm = nullptr) 
            : boundingBox(boundingBox), gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing), 
                minGlbCoordinate(boundingBox->getGlobalMin()), comm(comm) {
            
            const opp_point& minCoordinate = boundingBox->getGlobalMin();
            const opp_point& maxCoordinate = boundingBox->getGlobalMax();

            // this->globalGridDimensions.x = std::ceil((maxCoordinate.x - minCoordinate.x) * oneOverGridSpacing);
            // this->globalGridDimensions.y = std::ceil((maxCoordinate.y - minCoordinate.y) * oneOverGridSpacing);
            // this->globalGridDimensions.z = std::ceil((maxCoordinate.z - minCoordinate.z) * oneOverGridSpacing);   

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

            this->globalGridDimensions.x = ax;
            this->globalGridDimensions.y = ay;
            this->globalGridDimensions.z = az; 
        }

        //*******************************************************************************
        ~CellMapper() { 

#ifdef ENABLE_MPI
            // handle the shared memory window destroy routines
#else
            delete[] structMeshToCellMapping;
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
        inline int findClosestCellIndex(const opp_point& targetPosition) { 

            if (!boundingBox->isCoordinateInGlobalBoundingBox(targetPosition)) {
                // std::cerr << "isCoordinateInGlobalBoundingBox pos:" << targetPosition.x << "," << targetPosition.y << "," 
                //     << targetPosition.z << " is not within the bounding box" << std::endl;
                return MAX_CELL_INDEX;
            }

            size_t closestIndexMapping = getCellIndexMappingIndex(targetPosition);

            // Assume closestIndexMapping is within structMeshToCellMapping.size()
            return this->structMeshToCellMapping[closestIndexMapping];
        }

        //*******************************************************************************
        inline void generateStructMeshToCellIndexMapping(const opp_dat cellVolume_dat, const opp_dat cellDet_dat, 
            const opp_map cellConnectivity_map) { 
            
            size_t globalGridSize = (size_t)(globalGridDimensions.x * globalGridDimensions.y * globalGridDimensions.z);

            // if (OP_DEBUG)
                opp_printf("CellMapper", "generateStructMeshToCellIndexMapping global grid dimensions %d %d %d",
                    globalGridDimensions.x, globalGridDimensions.y, globalGridDimensions.z);

#ifdef ENABLE_MPI
            // create a shared memory window for everyone to write and then sync with all nodes
#else
            structMeshToCellMapping = new int[globalGridSize];
            for (size_t i = 0; i < globalGridSize; i++)
                structMeshToCellMapping[i] = MAX_CELL_INDEX;
#endif

            double lc[N_PER_C];
            int cell_set_size = cellVolume_dat->set->size;

            const opp_point& minCoordinate = boundingBox->getLocalMin();
            const opp_point& maxCoordinate = boundingBox->getLocalMax();

            double x = 0.0, y = 0.0, z = 0.0;

            for (int dz = 0; dz < globalGridDimensions.z; dz++) {
                
                z = minCoordinate.z + dz * this->gridSpacing;
                
                for (int dy = 0; dy < globalGridDimensions.y; dy++) {
                    
                    y = minCoordinate.y + dy * this->gridSpacing;
                    
                    for (int dx = 0; dx < globalGridDimensions.x; dx++) {
                        
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

                        if (cellIndex < cell_set_size) {
                            
                            // size_t index = getCellIndexMappingIndex(opp_point(x, y, z));
                            size_t index = (size_t)(dx + dy * globalGridDimensions.x + 
                                                    dz * globalGridDimensions.x * globalGridDimensions.y);
                            
                            this->structMeshToCellMapping[index] = cellIndex;
                        }  
                    }
                }
            }

            // printf("CellMapper %d\n", (int)globalGridSize);
            // for (size_t i = 0; i < globalGridSize; i++) {
            //     printf("%d ", this->structMeshToCellMapping[i]);

            //     if (i % 400 == 0) printf("\n");
            // }   
            // printf("\n");
        }

        // //*******************************************************************************
        // // This might have double precision issues
        // inline size_t getCellIndexMappingIndex(const opp_point& position) {

        //     int xIndex = static_cast<int>((position.x - this->minGlbCoordinate.x) * this->oneOverGridSpacing);
        //     int yIndex = static_cast<int>((position.y - this->minGlbCoordinate.y) * this->oneOverGridSpacing);
        //     int zIndex = static_cast<int>((position.z - this->minGlbCoordinate.z) * this->oneOverGridSpacing);

        //     return ((size_t)xIndex + (size_t)yIndex * globalGridDimensions.x + (size_t)zIndex * globalGridDimensions.x * globalGridDimensions.y);
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
            int xIndex = static_cast<int>(xDiff);
            int yIndex = static_cast<int>(yDiff);
            int zIndex = static_cast<int>(zDiff);

            // Calculate the cell index mapping index
            return ((size_t)(xIndex) + (size_t)(yIndex * globalGridDimensions.x) + 
                        (size_t)(zIndex * globalGridDimensions.x * globalGridDimensions.y));
        }

    private:
        const double gridSpacing = 0.0;
        const double oneOverGridSpacing = 0.0;
        const std::shared_ptr<BoundingBox> boundingBox;
        opp_ipoint globalGridDimensions;
        const opp_point& minGlbCoordinate;        
        int* structMeshToCellMapping = nullptr;         // This should contain mapping to global cell indices
        Comm* comm = nullptr;

        const double ONE_OVER_SIX = (1.0 / 6.0);
        const int N_PER_C = 4;
        const int DET_FIELDS = 4;
        const int NEIGHB_C = 4;
    };
};