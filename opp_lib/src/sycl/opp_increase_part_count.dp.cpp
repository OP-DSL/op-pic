
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

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <opp_sycl.h>

int* opp_saved_mesh_relation_d = nullptr;
size_t opp_saved_mesh_relation_size = 0;

//*************************************************************************************************
void opp_device_AssignMeshRelation(
    int *__restrict mesh_relation,
    const int *__restrict distribution,
    int start,
    int end,
    int inj_start,
    int inlet_size
    ,
    const sycl::nd_item<3> &item_ct1) 
{
    int tid = item_ct1.get_local_id(2) +
              item_ct1.get_group(2) * item_ct1.get_local_range(2);

    if (tid + start < end) 
    {    
        int n = tid + start;

        for (int i = 0; i < inlet_size; i++)
        {
            if (tid < distribution[i])
            {
                // assign inlet face index as the injected particle mesh relation
                mesh_relation[n + inj_start] = i; 
                break;
            } 
        }  
    }
}


//****************************************
void opp_inc_part_count_with_distribution(opp_set particles_set, 
    int num_particles_to_insert, opp_dat iface_dist, bool calc_new)
{
    if (OPP_DBG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d] %s", 
        num_particles_to_insert, (calc_new ? "NEW" : "COPY"));

    opp_profiler->start("IncPartCountWithDistribution");

    opp_dat mesh_rel_dat  = particles_set->mesh_relation_dat;

    // int nargs = 1;
    // opp_arg args[nargs];
    // args[0] = opp_arg_dat(mesh_rel_dat, OPP_READ);

    // int set_size = opp_mpi_halo_exchanges_grouped(particles_set, nargs, args, Device_CPU);

    // TODO : BUG What happens if the complete particle is dirty in device?

    opp_increase_particle_count(particles_set, num_particles_to_insert);

    int nargs1 = 2;
    opp_arg args1[nargs1];

    // if iface particle distribution is dirty in device, get it to the device
    args1[0] = opp_arg_dat(iface_dist, OPP_READ);
    args1[1] = opp_arg_dat(mesh_rel_dat, OPP_WRITE);

    int set_size = opp_mpi_halo_exchanges_grouped(particles_set, nargs1, args1, Device_GPU);
    opp_mpi_halo_wait_all(nargs1, args1);
    if (set_size > 0) 
    {
        int start     = 0;
        int end       = particles_set->diff;
        int inj_start = (particles_set->size - particles_set->diff);

        if (end - start > 0) 
        {
            if (calc_new) 
            {
                if (OPP_DBG) 
                    opp_printf("opp_inc_part_count_with_distribution", 
                        "Calculating all from new");

                int nthread = OPP_gpu_threads_per_block;
                int nblocks = (end - start - 1) / nthread + 1;

                /*
                DPCT1049:1: The work-group size passed to the SYCL kernel may
                exceed the limit. To get the device limit, query
                info::device::max_work_group_size. Adjust the work-group size if
                needed.
                */
                dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                    int *mesh_rel_dat_data_d_ct0 = (int *)mesh_rel_dat->data_d;
                    const int *iface_dist_data_d_ct1 =
                        (int *)iface_dist->data_d;
                    int iface_dist_set_size_ct5 = iface_dist->set->size;

                    cgh.parallel_for(
                        sycl::nd_range<3>(sycl::range<3>(1, 1, nblocks) *
                                              sycl::range<3>(1, 1, nthread),
                                          sycl::range<3>(1, 1, nthread)),
                        [=](sycl::nd_item<3> item_ct1) {
                            opp_device_AssignMeshRelation(
                                mesh_rel_dat_data_d_ct0, iface_dist_data_d_ct1,
                                start, end, inj_start, iface_dist_set_size_ct5,
                                item_ct1);
                        });
                });
            }
            else
            {
                size_t copy_size = (end - start) * sizeof(int);
                int* inj_mesh_relations = (int *)mesh_rel_dat->data_d + inj_start;

                if (opp_saved_mesh_relation_d == nullptr) 
                {
                    if (OPP_DBG) 
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "Allocating saved_mesh_relation_d with size [%zu]", copy_size);

                    opp_saved_mesh_relation_size = copy_size;
                    /*
                    DPCT1064:28: Migrated cudaMalloc call is used in a
                    macro/template definition and may not be valid for all
                    macro/template uses. Adjust the code.
                    */
                    cutilSafeCall(DPCT_CHECK_ERROR(
                        opp_saved_mesh_relation_d = (int *)sycl::malloc_device(
                            copy_size, dpct::get_in_order_queue())));
                    // opp_saved_mesh_relation_d = (int*)opp_host_malloc(copy_size);

                    // duplicate code below
                    int nthread = OPP_gpu_threads_per_block;
                    int nblocks = (end - start - 1) / nthread + 1;

                    /*
                    DPCT1049:2: The work-group size passed to the SYCL kernel
                    may exceed the limit. To get the device limit, query
                    info::device::max_work_group_size. Adjust the work-group
                    size if needed.
                    */
                    dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                        int *mesh_rel_dat_data_d_ct0 =
                            (int *)mesh_rel_dat->data_d;
                        const int *iface_dist_data_d_ct1 =
                            (int *)iface_dist->data_d;
                        int iface_dist_set_size_ct5 = iface_dist->set->size;

                        cgh.parallel_for(
                            sycl::nd_range<3>(sycl::range<3>(1, 1, nblocks) *
                                                  sycl::range<3>(1, 1, nthread),
                                              sycl::range<3>(1, 1, nthread)),
                            [=](sycl::nd_item<3> item_ct1) {
                                opp_device_AssignMeshRelation(
                                    mesh_rel_dat_data_d_ct0,
                                    iface_dist_data_d_ct1, start, end,
                                    inj_start, iface_dist_set_size_ct5,
                                    item_ct1);
                            });
                    });
                    // duplicate code above

                    cutilSafeCall(DPCT_CHECK_ERROR(
                        dpct::get_current_device().queues_wait_and_throw()));

                    cutilSafeCall(
                        DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(
                            opp_saved_mesh_relation_d, inj_mesh_relations,
                            copy_size)));
                }
                else
                {
                    if (OPP_DBG) 
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "Copying saved_mesh_relation_d with size [%zu]", copy_size);

                    if (opp_saved_mesh_relation_size != copy_size)
                    {
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "ERROR... saved_mesh_relation_size [%d] does not match with new copy size [%d]", 
                            opp_saved_mesh_relation_size, copy_size);
                    }

                    cutilSafeCall(
                        DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(
                            inj_mesh_relations, opp_saved_mesh_relation_d,
                            copy_size)));
                }
            }
        }
    }

    opp_set_dirtybit_grouped(nargs1, args1, Device_GPU);
    cutilSafeCall(
        DPCT_CHECK_ERROR(dpct::get_current_device().queues_wait_and_throw()));

    opp_profiler->end("IncPartCountWithDistribution");
}

