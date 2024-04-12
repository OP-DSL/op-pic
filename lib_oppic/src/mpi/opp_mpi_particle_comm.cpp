
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

#include <opp_mpi.h>

#define MPI_COUNT_EXCHANGE 0
#define MPI_TAG_PART_EX 1

#define PACK_AOS  false // PACK_AOS == false, is performing better

// this translate to std::map<oppic_set, std::map<local_cell_index, opp_particle_comm_data>>
std::map<oppic_set, std::map<int, opp_particle_comm_data>> opp_part_comm_neighbour_data; 

// this translate to std::map<particle_set, std::map<send_rank, std::vector<opp_particle_move_info>>>
std::map<int, std::map<int, std::vector<opp_particle_move_info>>> opp_part_move_indices;

std::vector<MPI_Request> send_req_count;
std::vector<MPI_Request> recv_req_count;

std::vector<int> opp_move_part_indices;

//*******************************************************************************
void opp_part_mark_move(oppic_set set, int particle_index, opp_particle_comm_data& comm_data)
{
    // if (OP_DEBUG) 
    //     opp_printf("opp_part_mark_move", "comm iter [%d] set [%s] | particle_index %d | send_rank %d | foreign_rank_index %d", 
    //         OPP_comm_iteration, set->name, particle_index, comm_data.cell_residing_rank, comm_data.local_index);

    auto& part_move_data_of_set = opp_part_move_indices[set->index];
    std::vector<opp_particle_move_info>& vec = part_move_data_of_set[comm_data.cell_residing_rank];

    opp_particle_move_info obj;
    obj.local_index = particle_index;
    obj.foreign_cell_index = comm_data.local_index;

    vec.push_back(obj);
}

//*******************************************************************************
void opp_part_pack(oppic_set set)
{
    if (OP_DEBUG) 
        opp_printf("opp_part_pack", "set [%s]", set->name);

    opp_profiler->start("Mv_Pack");

    opp_all_mpi_part_buffers* send_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    std::map<int, std::vector<opp_particle_move_info>>& set_move_indices_per_rank = opp_part_move_indices[set->index];

    // iterate over all the ranks to sent particles to
    for (auto& move_indices_per_rank : set_move_indices_per_rank)
    {
        int send_rank = move_indices_per_rank.first;
        std::vector<opp_particle_move_info>& move_indices_vector = move_indices_per_rank.second;

        opp_mpi_part_buffer& send_rank_buffer = send_buffers->buffers[send_rank];
        int64_t required_buffer_size = (move_indices_vector.size() * (int64_t)set->particle_size);

        // resize the export buffer if required
        if (send_rank_buffer.buf_export_index + required_buffer_size >= send_rank_buffer.buf_export_capacity)
        {
            if (send_rank_buffer.buf_export == nullptr)
            {
                send_rank_buffer.buf_export_capacity  = required_buffer_size;
                send_rank_buffer.buf_export_index     = 0;
                send_rank_buffer.buf_export           = (char *)opp_host_malloc(send_rank_buffer.buf_export_capacity);

                // opp_printf("opp_part_pack", "alloc buf_export capacity %" PRId64, send_rank_buffer.buf_export_capacity);
            }
            else 
            {
                // Assume that there are some particles left already, increase capacity beyond buf_export_index
                send_rank_buffer.buf_export_capacity  = send_rank_buffer.buf_export_index + required_buffer_size;
                send_rank_buffer.buf_export           = (char *)opp_host_realloc(send_rank_buffer.buf_export, 
                                                                        send_rank_buffer.buf_export_capacity);
                
                // opp_printf("opp_part_pack", "realloc buf_export capacity %" PRId64, send_rank_buffer.buf_export_capacity);
            }        
        }

#if PACK_AOS

        // iterate over all the particles to be sent to this send_rank
        for (auto& part : move_indices_vector)
        {
            int displacement = 0;

            // pack the particle data from dats into the export buffer 
            for (auto& dat : *(set->particle_dats))
            {
                if (dat->is_cell_index)
                {
                    // we need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &part.foreign_cell_index, dat->size);
                }
                else
                {
                    // copy the dat value to the send buffer
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &(dat->data[part.local_index * dat->size]), dat->size);
                }

                displacement += dat->size;
            } 

            send_rank_buffer.buf_export_index += set->particle_size;
            (send_buffers->export_counts)[send_rank] += 1;
        }

#else // PACK_SOA

        // Cannot use multiple packs before sending them, if opp_part_pack() is called multiple times with PACK_SOA, 
        // the communication data will get currupted
        int64_t displacement = 0;

        for (auto& dat : *(set->particle_dats))
        {
            int64_t dat_size = (int64_t)dat->size;

            if (dat->is_cell_index)
            {
                for (const auto& part : move_indices_vector)
                {
                    // we need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &part.foreign_cell_index, dat->size);
                 
                    displacement += dat_size;
                }

            }
            else
            {
                for (const auto& part : move_indices_vector)
                {
                    // copy the dat value to the send buffer
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &(dat->data[part.local_index * dat->size]), dat->size);
                    
                    displacement += dat_size;
                }                
            }
        }

        send_rank_buffer.buf_export_index = (int64_t)(set->particle_size * move_indices_vector.size()); // Not used
        (send_buffers->export_counts)[send_rank] = (int64_t)move_indices_vector.size();

#endif

        move_indices_vector.clear();
    }

    opp_profiler->end("Mv_Pack");

    if (OP_DEBUG) 
        opp_printf("opp_part_pack", "set [%s] END", set->name);
}

// //*******************************************************************************
// void opp_part_pack(oppic_set set, int index, int send_rank)
// {
//     // if (OP_DEBUG) 
//     //    opp_printf("opp_part_pack", "set [%s] | index %d | send_rank %d", set->name, index, send_rank);

//     opp_profiler->start("Mv_Pack");

//     opp_all_mpi_part_buffers* send_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

//     // check whether send_rank is a neighbour or not
//     if (OP_DEBUG)
//     {
//         if (std::find(send_buffers->neighbours.begin(), send_buffers->neighbours.end(), send_rank) 
//                         == send_buffers->neighbours.end())
//         {
//             opp_printf("opp_part_pack", "Error: send_rank %d is not a neighbour, cannot send index %d of set [%s]",
//                 send_rank, index, set->name);
//             opp_abort("opp_part_pack");
//         }
//     }

//     opp_mpi_part_buffer& send_rank_buffer = send_buffers->buffers[send_rank];

//     // resize the export buffer if required
//     if (send_rank_buffer.buf_export_index >= send_rank_buffer.buf_export_capacity)
//     {
//         if (send_rank_buffer.buf_export == nullptr)
//         {
//             send_rank_buffer.buf_export_capacity  = OPP_mpi_part_alloc_mult * set->particle_size;
//             send_rank_buffer.buf_export_index     = 0;
//             send_rank_buffer.buf_export           = (char *)opp_host_malloc(send_rank_buffer.buf_export_capacity);
//             //memset(send_rank_buffer.buf_export, 0, send_rank_buffer.buf_export_capacity); // not essential, can remove
//         }
//         else
//         {
//             send_rank_buffer.buf_export_capacity += OPP_mpi_part_alloc_mult * set->particle_size;
//             send_rank_buffer.buf_export           = (char *)opp_host_realloc(send_rank_buffer.buf_export, 
//                                                                     send_rank_buffer.buf_export_capacity);
//             // memset(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_capacity - 
//             //                                          OPP_mpi_part_alloc_mult * set->particle_size]),
//             //     0, OPP_mpi_part_alloc_mult * set->particle_size); // not essential, can remove
//         }

//         // opp_printf("opp_part_pack", "buf_export capacity %d", send_rank_buffer.buf_export_capacity);
//     }

//     std::vector<oppic_dat>& particle_dats = *(set->particle_dats);
//     int displacement = 0;

//     // pack the particle data from dats into the export buffer 
//     for (int i = 0; i < (int)particle_dats.size(); i++)
//     {
//         oppic_dat& dat = particle_dats[i];

//         memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
//                 &(dat->data[index * dat->size]), dat->size);

//         if (OP_DEBUG)
//         {
//             char* part_buffer = &send_rank_buffer.buf_export[send_rank_buffer.buf_export_index];
//             std::string log = "";
//             if (strcmp(dat->type, "double") == 0)
//             {
//                 double* d = (double*)(part_buffer + displacement);
//                 for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
//             }
//             else if (strcmp(dat->type, "int") == 0)
//             {
//                 int* d = (int*)(part_buffer + displacement);
//                 for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
//             }

//             // opp_printf("opp_part_pack", "%s from index %d -%s", dat->name, index, log.c_str());
//         }

//         displacement += dat->size;
//     } 

//     send_rank_buffer.buf_export_index += set->particle_size;
//     (send_buffers->export_counts)[send_rank] += 1;

//     opp_profiler->end("Mv_Pack");

//     // if (OP_DEBUG) 
//     //     opp_printf("opp_part_pack", "END send_rank %d exported count %d", send_rank, 
//     //         (send_buffers->export_counts)[send_rank]);
// }

//*******************************************************************************
void opp_part_unpack(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_unpack", "set [%s]", set->name);

    opp_profiler->start("Mv_Unpack");

    int64_t num_particles = 0;

    opp_all_mpi_part_buffers* recv_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;
    std::vector<int>& neighbours = recv_buffers->neighbours;

    // count the number of particles to be received from all ranks
    for (size_t i = 0; i < neighbours.size(); i++)
    {
        int neighbour_rank = neighbours[i];
        num_particles += (recv_buffers->import_counts)[neighbour_rank];
    }

    if (num_particles > 0)
    {
        if (!oppic_increase_particle_count_core(set, (int)num_particles)) // TODO : change this to int64_t
        {
            opp_printf("opp_part_unpack", "Error: Failed to increase particle count of particle set [%s]", set->name);
            opp_abort("opp_part_unpack");
        }

        int64_t new_part_index = (int64_t)(set->size - set->diff);

        for (int i = 0; i < (int)neighbours.size(); i++)
        {
            int recv_rank = neighbours[i];

            opp_mpi_part_buffer& receive_rank_buffer = recv_buffers->buffers[recv_rank];
            int64_t receive_count = recv_buffers->import_counts[recv_rank];

#if PACK_AOS
            std::vector<oppic_dat>& particle_dats = *(set->particle_dats);
            int64_t particle_size = set->particle_size;

            // unpack the received buffer from rank 'recv_rank' in to particle dats
            for (int part = 0; part < receive_count; part++)
            {
                char* part_buffer = &(receive_rank_buffer.buf_import[particle_size * part]);
                int displacement = 0;

                for (int i = 0; i < (int)particle_dats.size(); i++)
                {
                    oppic_dat& dat = particle_dats[i];

                    memcpy(dat->data + new_part_index * dat->size, part_buffer + displacement, dat->size);

                    displacement += dat->size;     
                }

                new_part_index++;
            }

#else // PACK_SOA

            int64_t displacement = 0;

            for (auto& dat : *(set->particle_dats))
            {
                memcpy(dat->data + new_part_index * dat->size, receive_rank_buffer.buf_import + displacement, 
                    dat->size * receive_count);
                
                displacement += dat->size * receive_count;            
            }

            new_part_index += receive_count;

#endif
        }
    }

    opp_profiler->end("Mv_Unpack");

    if (OP_DEBUG) opp_printf("opp_part_unpack", "END");
}

void opp_process_marked_particles(opp_set set) 
{
    for (int particle_index : opp_move_part_indices)
    {
        std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];
        const int map0idx = OPP_mesh_relation_data[particle_index];
        
        // Removed from the current rank
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;
        set->particle_remove_count += 1;

        auto it = set_part_com_data.find(map0idx);
        if (it == set_part_com_data.end())
        {
            opp_printf("opp_part_check_status", 
                "Error: cell %d cannot be found in opp_part_comm_neighbour_data map particle_index=%d", 
                    map0idx, particle_index);
            opp_abort("opp_part_check_status Error: cell cannot be found in opp_part_comm_neighbour_data map");
        }

        opp_particle_comm_data& comm_data = it->second;
        
        // change the cell_index to reflect the correct neighbour ranks local cell index before packing
        // OPP_mesh_relation_data[particle_index] = comm_data.local_index;
        // opp_part_pack(set, particle_index, comm_data.cell_residing_rank);

        opp_part_mark_move(set, particle_index, comm_data);
    }
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
bool opp_part_check_status(opp_move_var& m, int map0idx, oppic_set set, 
    int particle_index, int& remove_count, int thread) 
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        remove_count += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }
    else if (map0idx >= OPP_part_cells_set_size)
    {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        opp_move_part_indices.push_back(particle_index);
        return false;
    }

    // map0idx is an own cell and m.move_status == OPP_NEED_MOVE
    return true;
}

//*******************************************************************************
void opp_part_exchange(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_exchange", "set [%s] - particle size [%d]", set->name, set->particle_size);

    opp_profiler->start("Mv_Exchange");

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    std::vector<int>& neighbours = mpi_buffers->neighbours;
    int neighbour_count = neighbours.size();
    mpi_buffers->total_recv = 0;

    for (auto it = mpi_buffers->import_counts.begin(); it != mpi_buffers->import_counts.end(); it++)
        it->second = 0;

    mpi_buffers->recv_req.clear();
    mpi_buffers->send_req.clear();

    // std::vector<MPI_Request> send_req_count(neighbour_count);
    // std::vector<MPI_Request> recv_req_count(neighbour_count);
    send_req_count.clear();
    recv_req_count.clear();

    send_req_count.resize(neighbour_count);
    recv_req_count.resize(neighbour_count);

    opp_profiler->startMpiComm("", opp::OPP_Particle);

    // send/receive send_counts to/from only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];

        const int64_t& send_count = mpi_buffers->export_counts[neighbour_rank];
        MPI_Isend((void*)&send_count, 1, MPI_INT64_T, neighbour_rank, MPI_COUNT_EXCHANGE, OP_MPI_WORLD, &(send_req_count[i]));

        const int64_t& recv_count = mpi_buffers->import_counts[neighbour_rank];
        MPI_Irecv((void*)&recv_count, 1, MPI_INT64_T, neighbour_rank, MPI_COUNT_EXCHANGE, OP_MPI_WORLD, &(recv_req_count[i]));
    }

    double total_send_size = 0.0;

    // send the particle data only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];
        int64_t send_size = set->particle_size * mpi_buffers->export_counts[neighbour_rank];

        if (send_size <= 0)
        {
            if (OP_DEBUG) opp_printf("opp_part_exchange", "nothing to send to rank %d", neighbour_rank);
            continue;
        }
        else
        {   
            if (OP_DEBUG) opp_printf("opp_part_exchange", "sending %lld particle/s (size: %lld) to rank %d", 
                (send_size/ set->particle_size), send_size, neighbour_rank);
        }

        char* send_buffer = mpi_buffers->buffers[neighbour_rank].buf_export;

        // opp_printf("opp_part_exchange", "neighbour_rank: %d, send_size: %d", neighbour_rank, send_size);

        MPI_Request req;
        MPI_Isend(send_buffer, send_size, MPI_CHAR, neighbour_rank, MPI_TAG_PART_EX, OP_MPI_WORLD, &req);
        mpi_buffers->send_req.push_back(req);

        total_send_size += (send_size * 1.0f);
    }

    opp_profiler->addTransferSize("", opp::OPP_Particle, total_send_size, (size_t)(total_send_size / set->particle_size));
    opp_profiler->end("Mv_Exchange");

    std::string profName = std::string("Mv_WaitExCnt") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);
    // wait for the counts to receive only from neighbours
    MPI_Waitall(neighbour_count, &recv_req_count[0], MPI_STATUSES_IGNORE);
    opp_profiler->end(profName);

    opp_profiler->start("Mv_Exchange");
    // create/resize data structures and receive particle data from neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];
        int64_t recv_size = (int64_t)set->particle_size * mpi_buffers->import_counts[neighbour_rank];
        mpi_buffers->total_recv += mpi_buffers->import_counts[neighbour_rank];

        if (recv_size <= 0)
        {
            if (OP_DEBUG) opp_printf("opp_part_exchange", "nothing to receive from rank %d", neighbour_rank);
            continue;
        }


        opp_mpi_part_buffer& recv_buffer = mpi_buffers->buffers[neighbour_rank];

        if (recv_size >= recv_buffer.buf_import_capacity)
        {
            if (recv_buffer.buf_import == nullptr)
            {
                recv_buffer.buf_import_capacity  = recv_size;           
                recv_buffer.buf_import           = (char *)opp_host_malloc(recv_buffer.buf_import_capacity);          
            }
            else
            {
                recv_buffer.buf_import_capacity += recv_size;
                recv_buffer.buf_import           = (char *)opp_host_realloc(recv_buffer.buf_import, recv_buffer.buf_import_capacity);
            }
        }

        // opp_printf("opp_part_exchange", "neighbour_rank: %d, recv_size: %d, total_recv: %d, buff_import_capacity: %d", 
        //     neighbour_rank, recv_size, mpi_buffers->total_recv, recv_buffer.buf_import_capacity);

        MPI_Request req;
        MPI_Irecv(recv_buffer.buf_import, recv_size, MPI_CHAR, neighbour_rank, MPI_TAG_PART_EX, OP_MPI_WORLD, &req);
        mpi_buffers->recv_req.push_back(req);
    }

    // reset the export counts for another iteration
    for (auto it = mpi_buffers->export_counts.begin(); it != mpi_buffers->export_counts.end(); it++)
    {
        it->second = 0; // make the export count to zero for the next iteration
        mpi_buffers->buffers[it->first].buf_export_index = 0; // make the export index of that rank to zero for the next iteration
    }

    opp_profiler->endMpiComm("", opp::OPP_Particle);

    opp_profiler->end("Mv_Exchange");

    if (OP_DEBUG) opp_printf("opp_part_exchange", "END");
}

//*******************************************************************************
void cleanSendRecvBuffers(oppic_set set) 
{
    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    for (auto it = mpi_buffers->buffers.begin(); it != mpi_buffers->buffers.end(); it++)
    {
        opp_mpi_part_buffer& buffer = it->second;

        if (buffer.buf_import != nullptr)
        {
            buffer.buf_import_capacity = -1;  
            buffer.buf_import_index = 0;         
            opp_host_free(buffer.buf_import);
            buffer.buf_import = nullptr; 
        }

        if (buffer.buf_export != nullptr)
        {
            buffer.buf_export_capacity = -1;
            buffer.buf_export_index = 0;
            opp_host_free(buffer.buf_export);
            buffer.buf_export = nullptr;
        }
    }
}

//*******************************************************************************
void opp_part_wait_all(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_wait_all", "START");

    opp_profiler->start("Mv_WaitAll");

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    std::vector<MPI_Request>& send_req = mpi_buffers->send_req;
    std::vector<MPI_Request>& recv_req = mpi_buffers->recv_req;

    // std::vector<MPI_Status> recv_status(recv_req.size());
    // std::vector<MPI_Status> send_status(send_req.size());
    opp_profiler->startMpiComm("", opp::OPP_Particle);

    // wait till all the particles from all the ranks are received
    std::string profName = std::string("Mv_WaitExRecv") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);
    MPI_Waitall(send_req.size(), &(send_req[0]), MPI_STATUSES_IGNORE); // &(send_status[0])); //
    MPI_Waitall(recv_req.size(), &(recv_req[0]), MPI_STATUSES_IGNORE); // &(recv_status[0])); //
    opp_profiler->end(profName);

    opp_profiler->endMpiComm("", opp::OPP_Particle); // started at opp_part_exchange()

    send_req.clear();
    recv_req.clear();

    opp_profiler->end("Mv_WaitAll");

    if (OP_DEBUG) opp_printf("opp_part_wait_all", "END");
}

//*******************************************************************************
bool opp_part_check_all_done(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_check_all_done", "START");

    std::string profName = std::string("Mv_WaitDone") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;
    bool imported_parts = true;

    if (mpi_buffers->total_recv == 0)
    {
        imported_parts = false;
        set->diff = 0;
    }

    bool bool_ret = false;
    bool* buffer_recv = (bool *)opp_host_malloc(OPP_comm_size * sizeof(bool));

    opp_profiler->startMpiComm("", opp::OPP_Particle);

    // gather from all MPI ranks to see whether atleast one rank needs to communicate to another
    MPI_Allgather(&imported_parts, 1, MPI_C_BOOL, buffer_recv, 1, MPI_C_BOOL, OP_MPI_WORLD);

    opp_profiler->endMpiComm("", opp::OPP_Particle);

    std::string log = "";
    for (int i = 0; i < OPP_comm_size; i++)
    {
        bool_ret = (bool_ret || buffer_recv[i]);

        if (OP_DEBUG) 
            log += std::string(" R") + std::to_string(i) + (buffer_recv[i] ? "-T" : "-F");
    }

    if (OP_DEBUG) 
        opp_printf("opp_part_check_all_done", "iteration %d recv %lld - %s -%s", OPP_comm_iteration,
            mpi_buffers->total_recv, (bool_ret ? "ITER AGAIN" : "ALL DONE"), log.c_str());

    opp_host_free(buffer_recv);
    
    opp_profiler->end(profName);

    return !bool_ret;
}

//*******************************************************************************
void opp_part_comm_init()
{
    if (OP_DEBUG) opp_printf("opp_part_comm_init", "START");

    for (oppic_set& set : oppic_sets)
    {
        if (!set->is_particle) continue;

        opp_part_set_comm_init(set);
    }

    if (OP_DEBUG) opp_printf("opp_part_comm_init", "END");
}

//*******************************************************************************
void opp_part_set_comm_init(oppic_set set)
{
    // TODO : can use the same mappings for all particle sets with the same cells set, instead of communicating again

    if (OP_DEBUG) opp_printf("opp_part_set_comm_init", "set: %s", set->name);

    halo_list exp_exec_list = OP_export_exec_list[set->cells_set->index];
    halo_list imp_exec_list = OP_import_exec_list[set->cells_set->index];

    std::vector<MPI_Request> send_req;
    std::vector<MPI_Request> recv_req;    

    // send the local index of the export elements of the set to all neighbours
    for (int i = 0; i < exp_exec_list->ranks_size; i++) 
    {  
        int neighbour_rank = exp_exec_list->ranks[i];
        int* send_buffer = &(exp_exec_list->list[exp_exec_list->disps[i]]);
        int send_size = exp_exec_list->sizes[i];
    
        MPI_Request req;
        MPI_Isend(send_buffer, send_size, MPI_INT, neighbour_rank, OPP_rank, OP_MPI_WORLD, &req);
        send_req.push_back(req);  

        //print the per rank send buffers
        if (OP_DEBUG)
        {
            std::string log = "";
            for (int k = 0; k < send_size; k++)
                log += std::to_string(send_buffer[k]) + " ";        
            opp_printf("opp_part_set_comm_init", "%s SEND neighbour_rank %d send_size %d -> %s", 
                set->cells_set->name, neighbour_rank, send_size, log.c_str());
        } 
    }

    std::vector<std::vector<int>> recv_buffers(imp_exec_list->ranks_size);

    // receive the foreign ranks local index of the import elements of the set from all neighbours
    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = imp_exec_list->ranks[i];
        int recv_size = imp_exec_list->sizes[i];
        
        std::vector<int>& recv_buffer = recv_buffers[i];
        recv_buffer.resize(recv_size);
        
        MPI_Request req;
        MPI_Irecv(&recv_buffer[0], recv_size, MPI_INT, neighbour_rank, neighbour_rank, OP_MPI_WORLD, &req);
        recv_req.push_back(req);  
    }

    MPI_Waitall(recv_req.size(), &recv_req[0], MPI_STATUSES_IGNORE);
    
    // print the per rank received buffers
    if (OP_DEBUG)
    {
        for (int i = 0; i < (int)recv_buffers.size(); i++) 
        {
            // what I have (mappings) in import exec buffers, mappings before renumbering from that rank
            int* imp_buffer = &(imp_exec_list->list[imp_exec_list->disps[i]]); 

            std::string log = "";
            for (int k = 0; k < (int)recv_buffers[i].size(); k++)
                log += std::to_string((recv_buffers[i])[k]) + "|" + std::to_string(imp_buffer[k]) + " ";  

            opp_printf("opp_part_set_comm_init", "%s RECEIVE neighbour_rank %d recv_size %d (new|old) -> %s", 
                set->cells_set->name, imp_exec_list->ranks[i], recv_buffers[i].size(), log.c_str());
        }
    }

    if (true) // TODO : this might break existing OP2 functionality, check for issues
    { 
        for (int i = 0; i < (int)recv_buffers.size(); i++) 
        {
            int* imp_buffer = &(imp_exec_list->list[imp_exec_list->disps[i]]); 
            
            for (int k = 0; k < (int)recv_buffers[i].size(); k++)
                imp_buffer[k] = (recv_buffers[i])[k];  
        }
    }

    std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];

    // create mappings of neighbour ranks cell information for easy access during particle communication
    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = imp_exec_list->ranks[i];
        std::vector<int>& neighbour_rank_local_idxs = recv_buffers[i];

        for (int k = 0; k < imp_exec_list->sizes[i]; k++) 
        {       
            opp_particle_comm_data comm_data;
            comm_data.cell_residing_rank = neighbour_rank;
            comm_data.local_index = neighbour_rank_local_idxs[k];

            int local_index = (set->cells_set->size + imp_exec_list->disps[i] + k);

            set_part_com_data.insert({local_index, comm_data});

            // opp_printf("opp_part_comm_init", "set:[%s] li:[%d] nr:[%d] ni:[%d]", 
            //     set->cells_set->name, local_index, comm_data.cell_residing_rank, comm_data.local_index);
        }
    }  

    // create communication buffers in set
    opp_all_mpi_part_buffers* mpi_buffers = new opp_all_mpi_part_buffers();
    
    mpi_buffers->total_recv = 0;
    mpi_buffers->send_req.clear();
    mpi_buffers->recv_req.clear();

    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = imp_exec_list->ranks[i];

        mpi_buffers->neighbours.push_back(neighbour_rank);
        mpi_buffers->import_counts[neighbour_rank] = 0;

        opp_mpi_part_buffer& part_buffer = mpi_buffers->buffers[neighbour_rank];
        part_buffer.buf_import           = nullptr;
        part_buffer.buf_import_capacity  = -1;
        part_buffer.buf_import_index     = 0;
    }

    for (int i = 0; i < exp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = exp_exec_list->ranks[i];

        if (std::find(mpi_buffers->neighbours.begin(), mpi_buffers->neighbours.end(), neighbour_rank) 
                        == mpi_buffers->neighbours.end())
        {
            if (OP_DEBUG)
                opp_printf("opp_part_set_comm_init", "UNLIKELY, %d is on export list, but not on import list of set %s", 
                    neighbour_rank, set->name); 
            mpi_buffers->neighbours.push_back(neighbour_rank);
        }

        mpi_buffers->export_counts[neighbour_rank] = 0;

        opp_mpi_part_buffer& part_buffer = mpi_buffers->buffers[neighbour_rank];
        part_buffer.buf_export           = nullptr;
        part_buffer.buf_export_capacity  = -1;
        part_buffer.buf_export_index     = 0;
    }

    set->mpi_part_buffers = (void*)mpi_buffers;

    opp_profiler->reg("Mv_Pack");
    opp_profiler->reg("Mv_Exchange");
    opp_profiler->reg("Mv_Unpack");
    opp_profiler->reg("Mv_WaitAll");
    // opp_profiler->reg("Mv_WaitDone");
    // opp_profiler->reg("Mv_WaitExCnt");
    // opp_profiler->reg("Mv_WaitExRecv");

    std::string profName = "";
    for (int i = 0; i < 5; i++) {
        profName = std::string("Mv_WaitDone") + std::to_string(i);
        opp_profiler->reg(profName);
        profName = std::string("Mv_WaitExCnt") + std::to_string(i);
        opp_profiler->reg(profName);
        profName = std::string("Mv_WaitExRecv") + std::to_string(i);
        opp_profiler->reg(profName);
    }

    opp_part_move_indices.clear();

    if (OP_DEBUG) opp_printf("opp_part_set_comm_init", "set: %s END", set->name); 
}

//*******************************************************************************
void opp_part_comm_destroy()
{
    if (OP_DEBUG) opp_printf("opp_part_comm_destroy", "START"); 

    for (int i = 0; i < (int)oppic_sets.size(); i++)
    {
        oppic_set set = oppic_sets[i];

        if (!set->is_particle) continue;

        opp_all_mpi_part_buffers* all_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

        for (int rank : all_buffers->neighbours)
        {
            opp_mpi_part_buffer& buffer = all_buffers->buffers[rank];

            if (buffer.buf_export) 
                opp_host_free(buffer.buf_export);

            if (buffer.buf_import) 
                opp_host_free(buffer.buf_import);
        }

        delete all_buffers;
    }

    opp_part_comm_neighbour_data.clear();
    opp_part_move_indices.clear();

    if (OP_DEBUG) opp_printf("opp_part_comm_destroy", "END"); 
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
bool opp_part_checkForGlobalMove(opp_set set, const opp_point& point, const int partIndex, int& cellIdx) {
          
    size_t structCellIdx = cellMapper->findStructuredCellIndex(point);

    if (structCellIdx == MAX_CELL_INDEX) { // This happens when point is out of the unstructured mesh
        if (OP_DEBUG)
            opp_printf("opp_part_checkForGlobalMove", 
            "Remove %d [Struct cell index invalid - strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                partIndex, structCellIdx, point.x, point.y, point.z);
        cellIdx = MAX_CELL_INDEX;
        return true;
    }

    const int structCellRank = cellMapper->findClosestCellRank(structCellIdx);

    // Check whether the paticles need global moving, if yes start global moving process, 
    // if no, move to the closest local cell
    if (structCellRank != OPP_rank) {

        if (structCellRank == MAX_CELL_INDEX) {
            if (OP_DEBUG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [Rank invalid - strCellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, cellMapper->findClosestCellIndex(structCellIdx), structCellIdx, 
                    point.x, point.y, point.z);
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        const size_t globalCellIndex = cellMapper->findClosestCellIndex(structCellIdx);

        if (globalCellIndex == MAX_CELL_INDEX) {
            if (OP_DEBUG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [CellIdx invalid - strCellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, globalCellIndex, structCellIdx, point.x, point.y, point.z);
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // if the new rank is not the current rank, mark the particle to be sent via global comm
        globalMover->markParticleToMove(set, partIndex, structCellRank, globalCellIndex);

        // if (OP_DEBUG)
        //     opp_printf("opp_part_checkForGlobalMove", "Mark part %d [Move to rank %d gblCellIdx %d]", 
        //         partIndex, structCellRank, globalCellIndex);
        cellIdx = MAX_CELL_INDEX;
        return true;
    }
    else {
        
        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        cellIdx = cellMapper->findClosestCellIndex(structCellIdx);

        if (OP_DEBUG && (cellIdx < 0 || cellIdx >= set->cells_set->size)) {
            opp_printf("opp_part_checkForGlobalMove", 
                "Error... Particle %d assigned to current rank but invalid cell index %d [strCellIdx:%zu]", 
                    partIndex, cellIdx, structCellIdx);
            opp_abort("opp_part_checkForGlobalMove Error... Particle assigned to current rank but invalid cell index");
        }
    }
                
    return false;
}

//*******************************************************************************
//*******************************************************************************
using namespace opp;


ParticlePacker::ParticlePacker(std::map<int, std::map<int, std::vector<opp_particle_move_info>>>* partMoveData) 
    : partMoveData(partMoveData) {

}

//*******************************************************************************
ParticlePacker::~ParticlePacker() {

    this->partMoveData = nullptr;
    this->buffers.clear();
};

//*******************************************************************************
void ParticlePacker::pack(opp_set set) {

    if (OP_DEBUG) 
        opp_printf("ParticlePacker", "pack set [%s]", set->name);

    if (partMoveData == nullptr) {
        opp_abort(std::string("partMoveData is NULL in ParticlePacker"));
    }

    std::map<int, std::vector<char>>& buffersOfSet = this->buffers[set->index];
    // buffersOfSet.clear();
    for (auto& x : buffersOfSet) // try to keep the allocated vectors as it is, without deleting
        x.second.clear();

    if (partMoveData->at(set->index).size() == 0) {
        if (OP_DEBUG) 
            opp_printf("ParticlePacker", "Nothing to be sent for set [%s]", set->name);
        return;
    }     

    opp_profiler->start("GblMv_Pack");

    std::map<int, std::vector<opp_particle_move_info>>& moveIndicesPerSet = partMoveData->at(set->index);

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
                        &(dat->data[part.local_index * dat->size]), dat->size);
                    displacement += dat_size;
                }                
            }
        }

        if (OP_DEBUG)
            opp_printf("ParticlePacker", "Packed %zu parts to send to rank %d, displacement %d", 
                perRankBuffer.size(), send_rank, displacement);
    }

    opp_profiler->end("GblMv_Pack");
}

//*******************************************************************************
char* ParticlePacker::getBuffer(opp_set set, int sendRank) {

    auto itSet = this->buffers.find(set->index);

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
void ParticlePacker::unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                    int64_t totalParticlesToRecv, const std::vector<int64_t>& recvRankPartCounts) {

    if (OP_DEBUG) 
        opp_printf("ParticlePacker", "unpack set [%s]", set->name);

    opp_profiler->start("GblMv_Unpack");

    if (totalParticlesToRecv > 0)
    {
        std::vector<oppic_dat>& particle_dats = *(set->particle_dats);

        int64_t particle_size = set->particle_size;

        if (!oppic_increase_particle_count_core(set, (int)totalParticlesToRecv)) // TODO : make int to int64_t
        {
            opp_printf("ParticlePacker Unpack", "Error: Failed to increase particle count of particle set [%s]", 
                        set->name);
            opp_abort("ParticlePacker Unpack error");
        }

        int64_t newPartIndex = (int64_t)(set->size - set->diff);
        // int rankx = 0;

        for (const auto& x : particleRecvBuffers) {

            // int recvRank = x.first;
            const std::vector<char>& buffer = x.second;

            int64_t recvCount = ((int64_t)buffer.size() / particle_size) ; // recvRankPartCounts[rankx++];
            int64_t displacement = 0;

            for (auto& dat : particle_dats)
            {
                memcpy(dat->data + newPartIndex * dat->size, &(buffer[displacement]), 
                    dat->size * recvCount);
                
                displacement += dat->size * recvCount; //(dat->size * recvCount);            
            }

            newPartIndex += recvCount;
        }
    }

    opp_profiler->end("GblMv_Unpack");

    if (OP_DEBUG) 
        opp_printf("ParticlePacker", "Unpack END");
}


//*******************************************************************************
GlobalParticleMover::GlobalParticleMover(MPI_Comm comm) 
    : comm(comm) {

    CHECK(MPI_Win_allocate(sizeof(int), sizeof(int), MPI_INFO_NULL, this->comm,
                    &this->recv_win_data, &this->recv_win));
    
    globalPartMoveData.clear();

    packer = std::make_unique<ParticlePacker>(&globalPartMoveData);

    opp_profiler->reg("GblMv_WaitRanks");
    opp_profiler->reg("GblMv_Init");
    opp_profiler->reg("GblMv_Comm");
    opp_profiler->reg("GblMv_Finalize");
    opp_profiler->reg("GblMv_Pack");
    opp_profiler->reg("GblMv_Unpack");
    opp_profiler->reg("GblMv_WaitEx1");
    opp_profiler->reg("GblMv_WaitEx2");
    opp_profiler->reg("GblMv_WaitEx3");
    opp_profiler->reg("GblMv_WaitFin1");
    opp_profiler->reg("GblMv_WaitFin2");
}

//*******************************************************************************
GlobalParticleMover::~GlobalParticleMover() { 
    
    MPI_Barrier(this->comm);
    CHECK(MPI_Win_free(&this->recv_win));
}

//*******************************************************************************
inline void GlobalParticleMover::markParticleToMove(oppic_set set, int partIndex, int rankToBeMoved, int finalGlobalCellIndex) {
    
    // These validations should be already done
    // if (finalGlobalCellIndex == MAX_CELL_INDEX) {
    //     opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to rank %d but global index is invalid",
    //         partIndex, rankToBeMoved);
    //     return;
    // }

    // if (rankToBeMoved == MAX_CELL_INDEX) {
    //     opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to finalGlobalCellIndex %d but rank is invalid",
    //         partIndex, rankToBeMoved);
    //     return;
    // }

    auto& globalPartMoveDataOfSet = this->globalPartMoveData[set->index];
    std::vector<opp_particle_move_info>& vec = globalPartMoveDataOfSet[rankToBeMoved];

    opp_particle_move_info obj;
    obj.local_index = partIndex;
    obj.foreign_cell_index = finalGlobalCellIndex;
    vec.push_back(obj);
}

//*******************************************************************************
void GlobalParticleMover::initGlobalMove() {
    
    opp_profiler->start("GblMv_Init");

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

    opp_profiler->end("GblMv_Init");
}

//*******************************************************************************
void GlobalParticleMover::communicateParticleSendRecvRankCounts() {
    
    opp_profiler->start("GblMv_WaitRanks");

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

    opp_profiler->end("GblMv_WaitRanks");
}

//*******************************************************************************
void GlobalParticleMover::communicate(opp_set set) {
    
    opp_profiler->start("GblMv_Comm");

    std::map<int, std::vector<opp_particle_move_info>>& rankVsPartData = globalPartMoveData[set->index]; 
    this->numRemoteSendRanks = rankVsPartData.size();
    this->h_send_requests.resize(this->numRemoteSendRanks);
    this->h_send_ranks.resize(this->numRemoteSendRanks);
    this->h_send_rank_npart.resize(this->numRemoteSendRanks);

    int rankx = 0;
    for (auto& x : rankVsPartData) {

        if (x.first >= OPP_comm_size || x.first < 0) {
            opp_printf("GlobalParticleMover", "ERROR locking rank %d [size %zu] from rank %d", x.first, OPP_rank, x.second.size());
            this->numRemoteSendRanks -= 1;
            continue;
        }

        this->h_send_ranks[rankx] = x.first;
        this->h_send_rank_npart[rankx] = (int64_t)x.second.size();
        this->totalParticlesToSend += (int64_t)x.second.size();
        rankx++;
    }

    packer->pack(set);

    // (1) Gather the number of ranks that is going to send particles to the current rank ------------

    communicateParticleSendRecvRankCounts();

    opp_profiler->start("GblMv_WaitEx1");
    CHECK(MPI_Wait(&this->mpi_request, MPI_STATUS_IGNORE));
    opp_profiler->end("GblMv_WaitEx1");

    // At this point, the current rank knows how many particles to recv
    this->numRemoteRecvRanks = this->recv_win_data[0];

    this->h_recv_rank_npart.resize(this->numRemoteRecvRanks);
    this->h_recv_requests.resize(this->numRemoteRecvRanks);
    
    // (2) Exchange particle counts ----------------------------------------------------------------  

    // non-blocking recv of particle counts
    for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {

        CHECK(MPI_Irecv(&(this->h_recv_rank_npart[rankx]), 1, MPI_INT64_T, MPI_ANY_SOURCE, 
            42, this->comm, &(this->h_recv_requests[rankx])));
    }

    // non-blocking send of particle counts
    for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
        CHECK(MPI_Isend(&this->h_send_rank_npart[rankx], 1, MPI_INT64_T, this->h_send_ranks[rankx], 
            42, this->comm, &this->h_send_requests[rankx]));
    }

    this->h_recv_ranks.resize(this->numRemoteRecvRanks);
    this->h_recv_status.resize(this->numRemoteRecvRanks);

    opp_profiler->start("GblMv_WaitEx2");
    CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));
    opp_profiler->end("GblMv_WaitEx2");

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
        int64_t BytesToRecv = this->h_recv_rank_npart[rankx] * (int64_t)set->particle_size;
        partRecvBufferPerRank.resize(BytesToRecv);
    }

    if (OP_DEBUG)
        opp_printf("GlobalMove", "communicate %lld", this->totalParticlesToRecv);

    opp_profiler->start("GblMv_WaitEx3");
    CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));
    opp_profiler->end("GblMv_WaitEx3");

    // (3) send and receive the particles -----------------------------------------------------------
    
    this->h_send_requests.clear();
    this->h_recv_requests.clear();
    this->h_send_requests.resize(this->numRemoteSendRanks);            
    this->h_recv_requests.resize(this->numRemoteRecvRanks);
    
    for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {
        
        const int recvRank = this->h_recv_ranks[rankx];
        std::vector<char>& partRecvBufferPerRank = this->particleRecvBuffers[recvRank];
        int64_t recvSize = (this->h_recv_rank_npart[rankx] * set->particle_size);

        // opp_printf("Communicate", "Expected to recv %d bytes from rank %d, buffer size %zu", 
        //     sendSize, recvRank, partRecvBufferPerRank.size());

        CHECK(MPI_Irecv(&(partRecvBufferPerRank[0]), recvSize, MPI_CHAR, recvRank, 43, 
                this->comm, &this->h_recv_requests[rankx]));
    }

    // non-blocking send of particle data
    for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
        
        const int sendRank = this->h_send_ranks[rankx];
        char* sendBuffer = packer->getBuffer(set, sendRank);
        int64_t sendSize = (this->h_send_rank_npart[rankx] * set->particle_size);

        // opp_printf("Communicate", "Trying to send %d bytes to rank %d buffer %p", 
        //     sendSize, sendRank, sendBuffer);

        CHECK(MPI_Isend(sendBuffer, sendSize, MPI_CHAR, sendRank, 43, this->comm, 
                &this->h_send_requests[rankx]));
    }   

    // (4) Once sent, map could be cleared for the set, keep the allocations if possible -----------

    for (auto& x : globalPartMoveData[set->index])
        x.second.clear();       

    opp_profiler->end("GblMv_Comm");  
}

//*******************************************************************************
int64_t GlobalParticleMover::finalize(opp_set set) {
    
    opp_profiler->start("GblMv_Finalize");

    opp_profiler->start("GblMv_WaitFin1");
    CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));
    opp_profiler->end("GblMv_WaitFin1");

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

    opp_profiler->start("GblMv_WaitFin2");
    CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));
    opp_profiler->end("GblMv_WaitFin2");

    // Since packer->unpack might realloc dats, need to set the correct OPP_mesh_relation_data ptr
    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_profiler->end("GblMv_Finalize");

    return this->totalParticlesToRecv;
}