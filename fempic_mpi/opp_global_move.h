#pragma once

#include "opp_defs.h"
#include "opp_mpi.h"

namespace opp {

    //*******************************************************************************
    class ParticlePacker {

    public:
        //*******************************************************************************
        ParticlePacker(std::map<opp_set, std::map<int, std::vector<opp_particle_move_info>>>* partMoveData) 
            : partMoveData(partMoveData) {

        }

        //*******************************************************************************
        ~ParticlePacker() {

            this->partMoveData = nullptr;
            this->buffers.clear();
        };

        //*******************************************************************************
        inline void pack(opp_set set) {

            if (OP_DEBUG) 
                opp_printf("ParticlePacker", "pack set [%s]", set->name);

            if (partMoveData == nullptr) {
                opp_abort(std::string("partMoveData is NULL in ParticlePacker"));
            }

            std::map<int, std::vector<char>>& buffersOfSet = this->buffers[set];
            // buffersOfSet.clear();
            for (auto& x : buffersOfSet) // try to keep the allocated vectors as it is, without deleting
                x.second.clear();

            if (partMoveData->at(set).size() == 0) {
                if (OP_DEBUG) 
                    opp_printf("ParticlePacker", "Nothing to be sent for set [%s]", set->name);
                return;
            }     

            std::map<int, std::vector<opp_particle_move_info>>& moveIndicesPerSet = partMoveData->at(set);

            for (auto& moveIndicesPerRank : moveIndicesPerSet) {

                const int send_rank = moveIndicesPerRank.first;
                const std::vector<opp_particle_move_info>& moveIndicesVector = moveIndicesPerRank.second;

                size_t bytesPerRank = (size_t)set->particle_size * moveIndicesVector.size();
                
                std::vector<char>& perRankBuffer = buffersOfSet[send_rank];
                perRankBuffer.resize(bytesPerRank, 0);

                int displacement = 0;

                for (auto& dat : *(set->particle_dats)) {

                    int dat_size = dat->size;

                    if (dat->is_cell_index) {

                        for (const auto& part : moveIndicesVector) {

                            memcpy(&(perRankBuffer[displacement]), &part.foreign_cell_index, dat->size);
                            displacement += dat_size;
                        }
                    }
                    else {

                        for (const auto& part : moveIndicesVector) {

                            // copy the dat value to the send buffer
                            memcpy(&(perRankBuffer[displacement]), 
                                &(dat->data[part.local_particle_index * dat->size]), dat->size);
                            displacement += dat_size;
                        }                
                    }
                }

                if (OP_DEBUG)
                    opp_printf("ParticlePacker", "Packed %zu parts to send to rank %d, displacement %d", 
                        perRankBuffer.size(), send_rank, displacement);
            }
        }

        //*******************************************************************************
        inline char* getBuffer(opp_set set, int sendRank) {

            auto itSet = this->buffers.find(set);

            if (itSet == this->buffers.end()) {
                opp_abort(std::string("Set not found in ParticlePacker buffers"));
            } 

            auto itSetRank = itSet->second.find(sendRank);
            if (itSetRank == itSet->second.end()) {
                
                if (OP_DEBUG) 
                    opp_printf("ParticlePacker", "getBuffer set [%s] Rank [%d] does not have a buffer created", 
                        set->name, sendRank);
                return nullptr;
            } 

            if (OP_DEBUG)
                opp_printf("ParticlePacker", "getBuffer set [%s] Rank [%d] size %zu bytes", 
                        set->name, sendRank, itSetRank->second.size());

            return &(itSetRank->second[0]);
        }

        //*******************************************************************************
        inline void unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                            int totalParticlesToRecv, const std::vector<int>& recvRankPartCounts) {

            if (OP_DEBUG) 
                opp_printf("ParticlePacker", "unpack set [%s]", set->name);

            opp_profiler->start("Unpack");

            if (totalParticlesToRecv > 0)
            {
                std::vector<oppic_dat>& particle_dats = *(set->particle_dats);

                int particle_size = set->particle_size;

                if (!oppic_increase_particle_count_core(set, totalParticlesToRecv))
                {
                    opp_printf("ParticlePacker Unpack", "Error: Failed to increase particle count of particle set [%s]", 
                                set->name);
                    opp_abort("ParticlePacker Unpack error");
                }

                int newPartIndex = (set->size - set->diff);
                int rankx = 0;

                for (const auto& x : particleRecvBuffers) {

                    int recvRank = x.first;
                    const std::vector<char>& buffer = x.second;

                    int recvCount = ((int)buffer.size() / set->particle_size) ; // recvRankPartCounts[rankx++];
                    int displacement = 0;

                    for (auto& dat : particle_dats)
                    {
                        memcpy(dat->data + newPartIndex * dat->size, &(buffer[displacement]), 
                            dat->size * recvCount);
                        
                        displacement += dat->size * recvCount; //(dat->size * recvCount);            
                    }

                    newPartIndex += recvCount;
                }
            }

            opp_profiler->end("Unpack");

            if (OP_DEBUG) 
                opp_printf("ParticlePacker", "Unpack END");
        }

    private: 
        std::map<opp_set, std::map<int, std::vector<opp_particle_move_info>>>* partMoveData = nullptr;
        std::map<opp_set, std::map<int, std::vector<char>>> buffers;
    };


    //*******************************************************************************
    class GlobalParticleMover {

    public:
        //*******************************************************************************
        GlobalParticleMover(MPI_Comm comm) 
            : comm(comm) {

            CHECK(MPI_Win_allocate(sizeof(int), sizeof(int), MPI_INFO_NULL, this->comm,
                            &this->recv_win_data, &this->recv_win));
            
            globalPartMoveData.clear();

            packer = std::make_unique<ParticlePacker>(&globalPartMoveData);
        }

        //*******************************************************************************
        ~GlobalParticleMover() { 
            
            MPI_Barrier(this->comm);
            CHECK(MPI_Win_free(&this->recv_win));
        }

        //*******************************************************************************
        inline void markParticleToMove(oppic_set set, int partIndex, int rankToBeMoved, int finalGlobalCellIndex) {
            
            if (finalGlobalCellIndex == MAX_CELL_INDEX) {
                opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to rank %d but global index is invalid",
                    partIndex, rankToBeMoved);
                return;
            }

            if (rankToBeMoved == MAX_CELL_INDEX) {
                opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to finalGlobalCellIndex %d but rank is invalid",
                    partIndex, rankToBeMoved);
                return;
            }

            auto& globalPartMoveDataOfSet = this->globalPartMoveData[set];
            std::vector<opp_particle_move_info>& vec = globalPartMoveDataOfSet[rankToBeMoved];

            opp_particle_move_info obj;
            obj.local_particle_index = partIndex;
            obj.foreign_cell_index = finalGlobalCellIndex;
            vec.push_back(obj);
        }

        //*******************************************************************************
        inline void initGlobalMove() {
            
            this->h_send_ranks.clear();
            this->h_send_rank_npart.clear();
            this->h_recv_ranks.clear();
            this->h_recv_rank_npart.clear();
            this->h_send_requests.clear();
            this->h_recv_requests.clear();
            this->h_recv_status.clear();

            this->totalParticlesToSend = 0;
            this->totalParticlesToRecv = 0;
            this->numRemoteRecvRanks = 0;
            this->numRemoteSendRanks = 0;

            for (auto& x : particleRecvBuffers)
                x.second.clear();

            this->recv_win_data[0] = 0;
            CHECK(MPI_Ibarrier(this->comm, &this->mpi_request));
        }

        //*******************************************************************************
        inline void communicateParticleSendRecvRankCounts() {
           
            const int one[1] = {1};
            int recv[1];

            CHECK(MPI_Wait(&this->mpi_request, MPI_STATUS_IGNORE));

            for (auto& x : this->h_send_ranks) {
                
                const int rank = x;

                if (OP_DEBUG)
                    opp_printf("GlobalParticleMover", "locking rank %d", rank);

                CHECK(MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this->recv_win));
                CHECK(MPI_Get_accumulate(one, 1, MPI_INT, recv, 1, MPI_INT, rank, 0, 1,
                                            MPI_INT, MPI_SUM, this->recv_win));
                CHECK(MPI_Win_unlock(rank, this->recv_win));
            }

            CHECK(MPI_Ibarrier(this->comm, &this->mpi_request));   
        }

        //*******************************************************************************
        inline void communicate(opp_set set) {
            
            std::map<int, std::vector<opp_particle_move_info>>& rankVsPartData = globalPartMoveData[set]; 
            this->numRemoteSendRanks = rankVsPartData.size();
            this->h_send_requests.resize(this->numRemoteSendRanks);
            this->h_send_ranks.resize(this->numRemoteSendRanks);
            this->h_send_rank_npart.resize(this->numRemoteSendRanks);
     
            int rankx = 0;
            for (auto& x : rankVsPartData) {

                this->h_send_ranks[rankx] = x.first;
                this->h_send_rank_npart[rankx] = x.second.size();
                this->totalParticlesToSend += x.second.size();
                rankx++;
            }

            packer->pack(set);

            // (1) Gather the number of ranks that is going to send particles to the current rank ------------

            communicateParticleSendRecvRankCounts();

            CHECK(MPI_Wait(&this->mpi_request, MPI_STATUS_IGNORE));

            // At this point, the current rank knows how many particles to recv
            this->numRemoteRecvRanks = this->recv_win_data[0];

            this->h_recv_rank_npart.resize(this->numRemoteRecvRanks);
            this->h_recv_requests.resize(this->numRemoteRecvRanks);
            
            // (2) Exchange particle counts ----------------------------------------------------------------  

            // non-blocking recv of particle counts
            for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {

                CHECK(MPI_Irecv(&(this->h_recv_rank_npart[rankx]), 1, MPI_INT, MPI_ANY_SOURCE, 
                    42, this->comm, &(this->h_recv_requests[rankx])));
            }

            // non-blocking send of particle counts
            for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
                CHECK(MPI_Isend(&this->h_send_rank_npart[rankx], 1, MPI_INT, this->h_send_ranks[rankx], 
                    42, this->comm, &this->h_send_requests[rankx]));
            }

            this->h_recv_ranks.resize(this->numRemoteRecvRanks);
            this->h_recv_status.resize(this->numRemoteRecvRanks);

            CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));

            for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {

                const int remote_rank = this->h_recv_status[rankx].MPI_SOURCE;

                if (OP_DEBUG) { 

                    if (remote_rank < 0 && remote_rank < OPP_comm_size) 
                        opp_abort(std::string("Recv rank is invalid"));
                    if (this->h_recv_rank_npart[rankx] <= 0) 
                        opp_abort(std::string("A remote rank is trying to send 0 (or fewer) particles to this rank"));
                }

                this->h_recv_ranks[rankx] = remote_rank;
                this->totalParticlesToRecv += this->h_recv_rank_npart[rankx];

                std::vector<char>& partRecvBufferPerRank = particleRecvBuffers[remote_rank];
                size_t BytesToRecv = (size_t)this->h_recv_rank_npart[rankx] * set->particle_size;
                partRecvBufferPerRank.resize(BytesToRecv);
            }

            CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));

            // (3) send and receive the particles -----------------------------------------------------------
            
            this->h_send_requests.clear();
            this->h_recv_requests.clear();
            this->h_send_requests.resize(this->numRemoteSendRanks);            
            this->h_recv_requests.resize(this->numRemoteRecvRanks);
            
            for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {
                
                const int recvRank = this->h_recv_ranks[rankx];
                std::vector<char>& partRecvBufferPerRank = this->particleRecvBuffers[recvRank];
                int recvSize = (this->h_recv_rank_npart[rankx] * set->particle_size);

                // opp_printf("Communicate", "Expected to recv %d bytes from rank %d, buffer size %zu", 
                //     sendSize, recvRank, partRecvBufferPerRank.size());

                CHECK(MPI_Irecv(&(partRecvBufferPerRank[0]), recvSize, MPI_CHAR, recvRank, 43, 
                        this->comm, &this->h_recv_requests[rankx]));
            }

            // non-blocking send of particle data
            for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
                
                const int sendRank = this->h_send_ranks[rankx];
                char* sendBuffer = packer->getBuffer(set, sendRank);
                int sendSize = (this->h_send_rank_npart[rankx] * set->particle_size);

                // opp_printf("Communicate", "Trying to send %d bytes to rank %d buffer %p", 
                //     sendSize, sendRank, sendBuffer);

                CHECK(MPI_Isend(sendBuffer, sendSize, MPI_CHAR, sendRank, 43, this->comm, 
                        &this->h_send_requests[rankx]));
            }   

            // (4) Once sent, map could be cleared for the set, keep the allocations if possible -----------

            for (auto& x : globalPartMoveData[set])
                x.second.clear();       
        }

        //*******************************************************************************
        inline int finalize(opp_set set) {

            CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));

            // if (OP_DEBUG) 
            {
                // Check this rank recv'd the correct number of bytes from each remote
                for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {
                
                    MPI_Status *status = &this->h_recv_status[rankx];
                    int bytesRecvd = -1;
                    CHECK(MPI_Get_count(status, MPI_CHAR, &bytesRecvd));
                    int bytesExpected = (this->h_recv_rank_npart[rankx] * set->particle_size);

                    if (bytesRecvd != bytesExpected) {
                        
                        opp_printf("GlobalParticleMover", "recv'd incorrect num of bytes Expected %d Received %d",
                            bytesExpected, bytesRecvd);
                        opp_abort(std::string("recv'd incorrect number of bytes"));
                    }
                }
            }

            packer->unpack(set, this->particleRecvBuffers, this->totalParticlesToRecv, this->h_recv_rank_npart);

            // Note : The mesh relations received will be global cell indices and need to be converted to local

            CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));

            // Since packer->unpack might realloc dats, need to set the correct OPP_mesh_relation_data ptr
            OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

            return this->totalParticlesToRecv;
        }
    
    private:
        // this translate to std::map<particle_set, std::map<send_rank, std::vector<opp_particle_move_info>>>
        std::map<opp_set, std::map<int, std::vector<opp_particle_move_info>>> globalPartMoveData;
        MPI_Comm comm;

        int *recv_win_data;
        MPI_Win recv_win;
        MPI_Request mpi_request;

        std::vector<MPI_Request> h_send_requests;
        std::vector<MPI_Request> h_recv_requests;
        std::vector<MPI_Request> h_send_data_requests;
        std::vector<MPI_Request> h_recv_data_requests;
        std::vector<MPI_Status> h_recv_status;

        std::vector<int> h_send_ranks;
        std::vector<int> h_send_rank_npart;
        std::vector<int> h_recv_ranks;
        std::vector<int> h_recv_rank_npart;

        std::map<int, std::vector<char>> particleRecvBuffers; //std::map<int, std::pair<size_t, std::vector<char>>>

        int totalParticlesToSend = 0;
        int totalParticlesToRecv = 0;
        int numRemoteRecvRanks = 0;
        int numRemoteSendRanks = 0;

        std::unique_ptr<ParticlePacker> packer;
    };

} // namespace opp

