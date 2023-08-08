#pragma once

#include <mpi.h>
#include "opp_defs.h"

namespace opp {

    class Comm {

    public:    
        //*******************************************************************************
        Comm(MPI_Comm comm_parent) {
            this->comm_parent = comm_parent;

            int rank_parent;
            CHECK(MPI_Comm_rank(comm_parent, &rank_parent))
            CHECK(MPI_Comm_split_type(comm_parent, MPI_COMM_TYPE_SHARED, 0,
                                    MPI_INFO_NULL, &this->comm_intra))

            int rank_intra;
            CHECK(MPI_Comm_rank(this->comm_intra, &rank_intra))
            const int colour_intra = (rank_intra == 0) ? 1 : MPI_UNDEFINED;
            CHECK(MPI_Comm_split(comm_parent, colour_intra, 0, &this->comm_inter))

            CHECK(MPI_Comm_rank(this->comm_parent, &this->rank_parent))
            CHECK(MPI_Comm_rank(this->comm_intra, &this->rank_intra))
            CHECK(MPI_Comm_size(this->comm_parent, &this->size_parent))
            CHECK(MPI_Comm_size(this->comm_intra, &this->size_intra))
            
            if (comm_inter != MPI_COMM_NULL) {

                CHECK(MPI_Comm_rank(this->comm_inter, &this->rank_inter))
                CHECK(MPI_Comm_size(this->comm_inter, &this->size_inter))
            }

            opp_printf("Comm", "rank_parent %d|s=%d rank_intra %d|s=%d rank_inter %d|s=%d",
                this->rank_parent, this->size_parent, this->rank_intra, this->size_intra, 
                this->rank_inter, this->size_inter);
        };

        //*******************************************************************************
        ~Comm() {

            if ((this->comm_intra != MPI_COMM_NULL) && (this->comm_intra != MPI_COMM_WORLD)) {
                
                CHECK(MPI_Comm_free(&this->comm_intra))
                this->comm_intra = MPI_COMM_NULL;
            }

            if ((this->comm_inter != MPI_COMM_NULL) && (this->comm_inter != MPI_COMM_WORLD)) {

                CHECK(MPI_Comm_free(&this->comm_inter))
                this->comm_intra = MPI_COMM_NULL;
            }
        }

    public:
        /// Parent (i.e. global for the simulation) MPI communicator.
        MPI_Comm comm_parent;
        /// Communicator between one rank on each shared memory region.
        MPI_Comm comm_inter;
        /// Communicator between the ranks in a shared memory region.
        MPI_Comm comm_intra;
        /// MPI rank in the parent communicator.
        int rank_parent = -1;
        /// MPI rank in the inter shared memory region communicator.
        int rank_inter = -1;
        /// MPI rank within the shared memory region communicator.
        int rank_intra = -1;
        /// Size of the parent communicator.
        int size_parent = -1;
        /// Size of the inter shared memory communicator.
        int size_inter = -1;
        /// Size of the intra shared memory communicator.
        int size_intra = -1;
    };
};