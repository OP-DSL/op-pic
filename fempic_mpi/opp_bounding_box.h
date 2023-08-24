#pragma once

#include "opp_defs.h"
#include "oppic_lib_core.h"

#ifdef ENABLE_MPI
    #include "opp_comm.h"
#endif

#define TOLERENCE 1e-12

namespace opp {

    class BoundingBox {

    public:
        // For now, only 3D is implemented
        //*******************************************************************************
        BoundingBox(const opp_dat node_pos_dat, int dim, const std::shared_ptr<Comm> comm = nullptr) {
            
            if (dim != 3) {
                opp_abort(std::string("For now, only 3D BoundingBox is implemented"));
            }

            constexpr int DIM = 3;
            const double* node_pos_data = (const double*)node_pos_dat->data;
            const int node_pos_size = node_pos_dat->set->size + node_pos_dat->set->exec_size + 
                                        node_pos_dat->set->nonexec_size;;

            opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
            opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

            // make the bounding box even over the halo regions
            for (int i = 0; i < node_pos_size; i++) {
                minCoordinate.x = std::min(node_pos_data[i * DIM + 0], minCoordinate.x);
                minCoordinate.y = std::min(node_pos_data[i * DIM + 1], minCoordinate.y);
                minCoordinate.z = std::min(node_pos_data[i * DIM + 2], minCoordinate.z);
                maxCoordinate.x = std::max(node_pos_data[i * DIM + 0], maxCoordinate.x);
                maxCoordinate.y = std::max(node_pos_data[i * DIM + 1], maxCoordinate.y);
                maxCoordinate.z = std::max(node_pos_data[i * DIM + 2], maxCoordinate.z);
            }

            if (node_pos_size != 0) {
                minCoordinate.x -= TOLERENCE;
                minCoordinate.y -= TOLERENCE;
                minCoordinate.z -= TOLERENCE;
                maxCoordinate.x += TOLERENCE;
                maxCoordinate.y += TOLERENCE;
                maxCoordinate.z += TOLERENCE;
            }

            this->boundingBox[0] = minCoordinate;
            this->boundingBox[1] = maxCoordinate;

#ifdef ENABLE_MPI

            const double* localMin = reinterpret_cast<const double*>(&(this->boundingBox[0]));
            double* globalMin = reinterpret_cast<double*>(&(this->globalBoundingBox[0]));    
            MPI_Allreduce(localMin, globalMin, DIM, MPI_DOUBLE, MPI_MIN, comm->comm_parent);

            const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
            double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
            MPI_Allreduce(localMax, globalMax, DIM, MPI_DOUBLE, MPI_MAX, comm->comm_parent);

            // This is to avoid min and max corrdinates to not have MAX_REAL and MIN_REAL when current rank has no work
            if (node_pos_size == 0) {
                this->boundingBox[0].x = this->globalBoundingBox[0].x; 
                this->boundingBox[0].y = this->globalBoundingBox[0].y; 
                this->boundingBox[0].z = this->globalBoundingBox[0].z; 

                this->boundingBox[1].x = this->globalBoundingBox[0].x;
                this->boundingBox[1].y = this->globalBoundingBox[0].y;
                this->boundingBox[1].z = this->globalBoundingBox[0].z;
            }

            if (OPP_rank == OPP_ROOT)
                opp_printf("Global BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
                    this->globalBoundingBox[0].x, this->globalBoundingBox[0].y, this->globalBoundingBox[0].z, 
                    this->globalBoundingBox[1].x, this->globalBoundingBox[1].y, this->globalBoundingBox[1].z);
#else
            this->globalBoundingBox = this->boundingBox;
#endif

            if (OP_DEBUG)
                opp_printf("Local BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
                    this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
                    this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);
        }

        //*******************************************************************************
        ~BoundingBox() { };

        //*******************************************************************************
        const opp_point& getLocalMin() const {

            return this->boundingBox[0];
        }

        //*******************************************************************************
        const opp_point& getLocalMax() const {

            return this->boundingBox[1];
        }

        //*******************************************************************************
        const opp_point& getGlobalMin() const {

            return this->globalBoundingBox[0];
        }

        //*******************************************************************************
        const opp_point& getGlobalMax() const {
            
            return this->globalBoundingBox[1];
        }

        //*******************************************************************************
        bool isCoordinateInBoundingBox(const opp_point& point) { 

            if (this->boundingBox[0].x > point.x || this->boundingBox[1].x < point.x) 
                return false;
            else if (this->boundingBox[0].y > point.y || this->boundingBox[1].y < point.y) 
                return false;
            else if (this->boundingBox[0].z > point.z || this->boundingBox[1].z < point.z) 
                return false;
            
            return true;
        }

        //*******************************************************************************
        bool isCoordinateInGlobalBoundingBox(const opp_point& point) { 

            if (this->globalBoundingBox[0].x > point.x || this->globalBoundingBox[1].x < point.x) 
                return false;
            else if (this->globalBoundingBox[0].y > point.y || this->globalBoundingBox[1].y < point.y) 
                return false;
            else if (this->globalBoundingBox[0].z > point.z || this->globalBoundingBox[1].z < point.z) 
                return false;
            
            return true;
        }

    private:
        std::array<opp_point,2> boundingBox; // index 0 is min, index 1 is max
        std::array<opp_point,2> globalBoundingBox; // index 0 is min, index 1 is max
    };
};