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

// #pragma once

#include <opp_cuda.h>

//*******************************************************************************
void opp_part_pack_cuda(oppic_set set);
void opp_part_unpack_cuda(oppic_set set);
void particle_sort_cuda(oppic_set set, bool hole_filling);

std::vector<char> OPP_need_remove_flags;
char *OPP_need_remove_flags_d = nullptr;
int OPP_need_remove_flags_size = 0;

thrust::host_vector<int> OPP_thrust_move_indices_h;
thrust::device_vector<int> OPP_thrust_move_indices_d;
int *OPP_move_indices_d = nullptr;
int *OPP_move_count_d = nullptr;
int OPP_move_count_h = 0;

// int OPP_move_indices_capacity = 0;

struct CopyMaxCellIndexFunctor {
    int* B;
    CopyMaxCellIndexFunctor(int* B) : B(B) {}

    __host__ __device__
    void operator()(int index) const {
        B[index] = MAX_CELL_INDEX;
    }
};

__global__ void setArrayToMaxCID(int* array, int size) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < size) {
        array[tid] = MAX_CELL_INDEX;
    }
}

//*******************************************************************************
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args)
{ 

    oppic_init_particle_move_core(set);

    cutilSafeCall(cudaMemcpy(set->particle_remove_count_d, &(set->particle_remove_count), sizeof(int), 
                    cudaMemcpyHostToDevice));

    // if (set->size >= OPP_need_remove_flags_size)
    // {
    //     if (OPP_need_remove_flags_d != NULL) 
    //         cutilSafeCall(cudaFree(OPP_need_remove_flags_d));
        
    //     OPP_need_remove_flags_size = set->size;
    //     cutilSafeCall(cudaMalloc(&OPP_need_remove_flags_d, OPP_need_remove_flags_size * sizeof(char)));
    // }

    // cutilSafeCall(cudaMemset(OPP_need_remove_flags_d, 0, OPP_need_remove_flags_size));

    if (set->size > (int)OPP_thrust_move_indices_d.size())
    {     
        OPP_thrust_move_indices_d.resize(set->size);
        OPP_thrust_move_indices_h.resize(set->size);

        OPP_move_indices_d = (int*)thrust::raw_pointer_cast(OPP_thrust_move_indices_d.data());

        if (OPP_move_count_d == nullptr) {
            cutilSafeCall(cudaMalloc(&OPP_move_count_d, sizeof(int)));
        }
    }

    OPP_move_count_h = 0;
    cutilSafeCall(cudaMemcpy(OPP_move_count_d, &OPP_move_count_h, sizeof(int), cudaMemcpyHostToDevice));

    // TODO : redundant - remove below
    // init the array to MAX_CELL_INDEX
    // const int blocksPerGrid = (size - 1) / OPP_gpu_threads_per_block + 1;
    // setArrayToMaxCID<<<blocksPerGrid, OPP_gpu_threads_per_block>>>(OPP_move_indices_d, set->size);

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }
    else
    {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++)
        {
            if (args[i].argtype == OP_ARG_DAT && args[i].dat->set->is_particle)
            {
                args[i].data = args[i].dat->data;
                args[i].data_d = args[i].dat->data_d;
                if (OP_DEBUG) opp_printf("SSSS", "dat %s %p", args[i].dat->name, args[i].dat->data_d);
            }
        }
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
    OPP_mesh_relation_data_d = ((int *)set->mesh_relation_dat->data_d); 
}

//*******************************************************************************
// 1. keep track of whether atleast one need mpi comm [WRONG] I might not require to send to others, but some one might try to send to me 
// 2. if yes, download dats from device and use opp_part_exchange (in this particle will get copied so we can override the space)
// 3. do the below, oppic_particle_sort or particle_sort_cuda
// 4. check whether all mpi ranks are done and if yes, return false
// 5. if no, wait for all to complete and copy only the received particles to the device buffer 
//        (may need to write another path and to expand device array sizes)
// 6. set new start and ends and return true
bool opp_finalize_particle_move(oppic_set set)
{ 
    opp_profiler->start("Mv_Finalize");

    cutilSafeCall(cudaDeviceSynchronize());

    OPP_move_count_h = 0;
    cutilSafeCall(cudaMemcpy(&OPP_move_count_h, OPP_move_count_d, sizeof(int), 
        cudaMemcpyDeviceToHost));

    cutilSafeCall(cudaMemcpy(&(set->particle_remove_count), set->particle_remove_count_d, 
                    sizeof(int), cudaMemcpyDeviceToHost));

    const int particle_remove_count = set->particle_remove_count;

    if (OP_DEBUG) 
        opp_printf("oppic_finalize_particle_move", "set [%s][%d] remove_count [%d] move count [%d]", 
            set->name, set->size, particle_remove_count, OPP_move_count_h);

#ifdef USE_MPI
    // At this stage, particles of device is clean

    // download only the required particles to send and pack them in rank based mpi buffers
    opp_profiler->start("Mv_F_pack");
    opp_part_pack_cuda(set);
    opp_profiler->end("Mv_F_pack");

    // send the counts and send the particles  
    opp_profiler->start("Mv_F_ex");   
    opp_part_exchange(set); 
    opp_profiler->end("Mv_F_ex");

#endif

    opp_profiler->start("Mv_F_fill");
    if (particle_remove_count > 0)
    {
        set->size -= particle_remove_count;

        if (OP_auto_sort == 1)
        {
            if (OP_DEBUG) 
                opp_printf("oppic_finalize_particle_move", "auto sorting particle set [%s]", 
                set->name);
            oppic_particle_sort(set);
        }
        else
        {
            particle_sort_cuda(set, true); // Does only hole filling
        }
    }
    opp_profiler->end("Mv_F_fill");

#ifdef USE_MPI
    opp_profiler->start("Mv_F_check");
    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");
        opp_profiler->end("Mv_F_check");
        return false; // all mpi ranks do not have anything to communicate to any rank
    }
    opp_profiler->end("Mv_F_check");

    opp_profiler->start("Mv_F_wait");
    opp_part_wait_all(set); // wait till all the particles are communicated
    opp_profiler->end("Mv_F_wait");

    if (OP_DEBUG)
        opp_printf("opp_finalize_particle_move", "set [%s] size prior unpack %d", set->name, set->size);

    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_profiler->start("Mv_F_unpack");
    opp_part_unpack_cuda(set);    
    opp_profiler->end("Mv_F_unpack");

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");

    return true;
#else
    return false;
#endif
}

// Cannot use multiple packs before sending them, if opp_part_pack() is called multiple times with PACK_SOA, 
// the communication data may get currupted
//*******************************************************************************
void opp_part_pack_cuda(opp_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_pack_cuda", "start");

#ifdef USE_MPI
    opp_profiler->start("Mv_Pack");

//opp_profiler->start("Mv_Pack_resize");
    // OPP_need_remove_flags.resize(OPP_need_remove_flags_size, 0);
//opp_profiler->end("Mv_Pack_resize");

//opp_profiler->start("Mv_Pack_d_to_h");
    // cutilSafeCall(cudaMemcpy((char*)&(OPP_need_remove_flags[0]), OPP_need_remove_flags_d, 
    //                 OPP_need_remove_flags_size, cudaMemcpyDeviceToHost));
//opp_profiler->end("Mv_Pack_d_to_h");

//opp_profiler->start("Mv_Pack_gather");
    // gather the particle indices to be sent
    // thrust::host_vector<int> send_indices_hv;
    // for (size_t index = 0; index < OPP_need_remove_flags.size(); index++) {
    //     if (OPP_need_remove_flags[index] == 1) {
    //         send_indices_hv.push_back(index);
    //     }
    // }
    // auto findIfPredicate = [](char value) { return value == 1; };
    // auto findIfBegin = std::find_if(OPP_need_remove_flags.begin(), OPP_need_remove_flags.end(), findIfPredicate);
    // auto findIfEnd = OPP_need_remove_flags.end();
    // while (findIfBegin != findIfEnd) 
    // {
    //     send_indices_hv.push_back(std::distance(OPP_need_remove_flags.begin(), findIfBegin));
    //     findIfBegin = std::find_if(std::next(findIfBegin), findIfEnd, findIfPredicate);
    // }

//opp_profiler->end("Mv_Pack_gather");

    // int part_send_count = (int)send_indices_hv.size();
    if (OPP_move_count_h <= 0) 
    {
        opp_profiler->end("Mv_Pack");
        return;
    }

    thrust::sort(OPP_thrust_move_indices_d.begin(), 
        OPP_thrust_move_indices_d.begin() + OPP_move_count_h);

    // thrust::device_vector<int> send_indices_dv(send_indices_hv);
    thrust::copy(OPP_thrust_move_indices_d.begin(), 
        OPP_thrust_move_indices_d.begin() + OPP_move_count_h, OPP_thrust_move_indices_h.begin());
    
    // std::string loggg = std::to_string(OPP_move_count_h) + " | ";
    // for (int k = 0; k < OPP_move_count_h; k++)
    //     loggg += std::to_string(OPP_thrust_move_indices_h[k]) + std::string(" ");
    // opp_printf("opp_part_pack_cuda", "Part indices : %s", loggg.c_str());

    // copy the cell indices of the particles to be sent
    thrust::device_vector<int> send_part_cell_idx_dv(OPP_move_count_h);
    copy_according_to_index(set->mesh_relation_dat->thrust_int, &send_part_cell_idx_dv, 
        OPP_thrust_move_indices_d, -1, -1, OPP_move_count_h, 1);
    thrust::host_vector<int> send_part_cell_idx_hv(send_part_cell_idx_dv);

    // enrich the particles to communicate with the correct external cell index and mpi rank
    std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];
    for (int index = 0; index < OPP_move_count_h; index++)
    {
        int particle_index = OPP_thrust_move_indices_h[index];
        int map0idx = send_part_cell_idx_hv[index];

        auto it = set_part_com_data.find(map0idx);
        if (it == set_part_com_data.end())
        {
            opp_printf("opp_part_pack_cuda", 
                "Error: cell %d cannot be found in opp_part_comm_neighbour_data map", map0idx);
            return; // unlikely, need opp_abort() instead!
        }

        opp_part_mark_move(set, index, it->second);
    }

    std::map<int, std::vector<char>> move_dat_data_map;

    // download the particles to send
    {
        thrust::device_vector<int> temp_int_dv;
        thrust::device_vector<double> temp_real_dv;

        for (auto& dat : *(set->particle_dats)) 
        {
            size_t bytes_to_copy = (OPP_move_count_h * dat->size);
            
            auto& move_dat_data = move_dat_data_map[dat->index];
            move_dat_data.resize(bytes_to_copy);

            if (strcmp(dat->type, "double") == 0)
            {
                temp_real_dv.resize(OPP_move_count_h * dat->dim);
                copy_according_to_index(dat->thrust_real, &temp_real_dv, OPP_thrust_move_indices_d, 
                    dat->set->set_capacity, OPP_move_count_h, OPP_move_count_h, dat->dim);
                
                cudaMemcpy(move_dat_data.data(), thrust::raw_pointer_cast(&temp_real_dv[0]), 
                    bytes_to_copy, cudaMemcpyDeviceToHost);
            }
            else if (strcmp(dat->type, "int") == 0)
            {
                temp_int_dv.resize(OPP_move_count_h * dat->dim);
                copy_according_to_index(dat->thrust_int, &temp_int_dv, OPP_thrust_move_indices_d, 
                    dat->set->set_capacity, OPP_move_count_h, OPP_move_count_h, dat->dim);
                
                cudaMemcpy(move_dat_data.data(), thrust::raw_pointer_cast(&temp_int_dv[0]), 
                    bytes_to_copy, cudaMemcpyDeviceToHost);
            }
        }      
    }

    opp_all_mpi_part_buffers* send_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    // increase the sizes of MPI buffers
    for (auto& move_indices_per_rank : opp_part_move_indices[set->index])
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
                send_rank_buffer.buf_export_capacity  = OPP_mpi_part_alloc_mult * required_buffer_size;
                send_rank_buffer.buf_export_index     = 0;
                send_rank_buffer.buf_export           = (char *)malloc(send_rank_buffer.buf_export_capacity);

                // opp_printf("opp_part_pack", "alloc buf_export capacity %d", send_rank_buffer.buf_export_capacity);
            }
            else 
            {
                // Assume that there are some particles left already, increase capacity beyond buf_export_index
                send_rank_buffer.buf_export_capacity  = send_rank_buffer.buf_export_index + 
                                                            OPP_mpi_part_alloc_mult * required_buffer_size;
                send_rank_buffer.buf_export           = (char *)realloc(send_rank_buffer.buf_export, 
                                                            send_rank_buffer.buf_export_capacity);
                
                // opp_printf("opp_part_pack", "realloc buf_export capacity %d", send_rank_buffer.buf_export_capacity);
            }        
        }
    }

    // iterate over all the ranks and pack to mpi buffers using SOA
    for (auto& move_indices_per_rank : opp_part_move_indices[set->index])
    {
        int send_rank = move_indices_per_rank.first;
        std::vector<opp_particle_move_info>& move_indices_vector = move_indices_per_rank.second;
        size_t per_rank_move_count = move_indices_vector.size();

        opp_mpi_part_buffer& send_rank_buffer = send_buffers->buffers[send_rank];

        int64_t displacement = 0;
        for (auto& dat : *(set->particle_dats))
        {
            auto& move_dat_data = move_dat_data_map[dat->index];

            int64_t dat_size = (int64_t)dat->size;
            int64_t element_size = (int64_t)(dat->size / dat->dim);

            if (dat->is_cell_index)
            {
                for (const auto& move_info : move_indices_vector)
                {
                    // we need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &move_info.foreign_cell_index, dat->size);
                 
                    displacement += dat_size;
                }
            }
            else
            {
                for (int d = 0; d < dat->dim; d++) 
                {
                    for (const auto& move_info : move_indices_vector)
                    {
                        // copy the multi dimensional dat value to the send buffer
                        memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                            &(move_dat_data[(d * OPP_move_count_h + move_info.local_index) * element_size]), 
                            element_size);
                        
                        displacement += element_size;
                    }
                }                
            }
        }

        send_rank_buffer.buf_export_index = (int64_t)(set->particle_size * move_indices_vector.size()); // Not used
        (send_buffers->export_counts)[send_rank] = (int64_t)move_indices_vector.size();

        move_indices_vector.clear();
    }

    // This particle is already packed, hence needs to be removed from the current rank
    CopyMaxCellIndexFunctor copyMaxCellIndexFunctor((int*)set->mesh_relation_dat->data_d);
    thrust::for_each(OPP_thrust_move_indices_d.begin(), OPP_thrust_move_indices_d.begin() + OPP_move_count_h, 
        copyMaxCellIndexFunctor);

    opp_profiler->end("Mv_Pack");
#endif

    if (OP_DEBUG) opp_printf("opp_part_pack_cuda", "end");
}

//*******************************************************************************
void opp_part_unpack_cuda(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_part_unpack_cuda", "set [%s]", set->name);

#ifdef USE_MPI
    opp_profiler->start("Mv_Unpack");

    std::vector<oppic_dat>& particle_dats = *(set->particle_dats);
    int64_t num_new_particles = 0;

    opp_all_mpi_part_buffers* recv_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;
    std::vector<int>& neighbours = recv_buffers->neighbours;

    // count the number of particles to be received from all ranks
    for (size_t i = 0; i < neighbours.size(); i++)
    {
        int neighbour_rank = neighbours[i];
        num_new_particles += (recv_buffers->import_counts)[neighbour_rank];
    }

    if (num_new_particles > 0)
    {
        int64_t new_part_index = (int64_t)(set->size);

        oppic_increase_particle_count(set, (int)num_new_particles);

        // create a continuous memory in host to copy to device
        std::vector<char*> temp_data_vec;
        for (size_t d = 0; d < particle_dats.size(); d++)
        {
            opp_dat dat = particle_dats[d];
            char *temp_data = (char *)malloc(dat->size * num_new_particles * sizeof(char));
            temp_data_vec.push_back(temp_data);
        }

        int current_recv_count = 0;
        for (int i = 0; i < (int)neighbours.size(); i++)
        {
            int recv_rank = neighbours[i];

            opp_mpi_part_buffer& receive_rank_buffer = recv_buffers->buffers[recv_rank];

            int64_t receive_count = recv_buffers->import_counts[recv_rank];
            int64_t displacement = 0;

            for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++)
            {
                opp_dat dat = particle_dats[dat_idx];
                char *temp_data = temp_data_vec[dat_idx];

                int element_size = dat->size / dat->dim;

                for (int i = 0; i < dat->dim; i++) 
                {
                    memcpy(&(temp_data[element_size * (i * num_new_particles + current_recv_count)]), 
                        &(receive_rank_buffer.buf_import[displacement]), element_size * receive_count);

                    displacement += element_size * receive_count; 
                }                
            }

            // for (size_t d = 0; d < particle_dats.size(); d++)
            // {
            //     opp_dat dat = particle_dats[d];
            //     char *temp_data = temp_data_vec[d];

            //     int element_size = dat->size / dat->dim;

            //     for (int i = 0; i < dat->dim; i++) 
            //     {
            //         for (int64_t j = 0; j < receive_count; j++) 
            //         {
            //             int64_t tmp_index = element_size * i * num_new_particles + element_size * (j + current_recv_count);
            //             int64_t recv_index = displacement + dat->size * j + element_size * i;

            //             for (int c = 0; c < element_size; c++) 
            //             {
            //                 temp_data[tmp_index + c] = receive_rank_buffer.buf_import[recv_index + c];
            //             }
            //         }
            //     }

            //     displacement += dat->size * receive_count;               
            // }

            // current_recv_count += receive_count;

            current_recv_count += receive_count;
        }

        // copy to device
        for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++)
        {
            opp_dat dat = particle_dats[dat_idx];
            char *temp_data = temp_data_vec[dat_idx];

            size_t bytes_to_copy_per_dim = num_new_particles * dat->size / dat->dim;
            int element_size = dat->size / dat->dim;

            for (int64_t d = 0; d < dat->dim; d++) 
            {     
                size_t data_d_offset = (new_part_index + d * set->set_capacity) * element_size;
                size_t data_h_offset = d * num_new_particles * element_size;
    
                cutilSafeCall(cudaMemcpy((dat->data_d + data_d_offset), (temp_data + data_h_offset), 
                                    bytes_to_copy_per_dim, cudaMemcpyHostToDevice));       
            }
        }

        cutilSafeCall(cudaDeviceSynchronize());

        for (auto& char_ptr : temp_data_vec)
            free(char_ptr);
        temp_data_vec.clear();
    }

    opp_profiler->end("Mv_Unpack");
#endif

    if (OP_DEBUG) opp_printf("opp_part_unpack_cuda", "END");    
}

//*******************************************************************************