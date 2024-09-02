//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;

OPP_INT* opp_k1_dat0_stride_s = nullptr;
OPP_INT* opp_k1_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__update_pos_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1 // p_pos | OPP_RW
)
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) 
        opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    
    opp_set_stride(opp_k1_dat0_stride_s, opp_k1_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat1_stride_s, opp_k1_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        num_blocks = (end - start - 1) / block_size + 1;

        {
            opp_queue->submit([&](sycl::handler &cgh) {

                OPP_REAL* arg0_data_d = (OPP_REAL *)args[0].data_d;     // p_vel
                OPP_REAL* arg1_data_d = (OPP_REAL *)args[1].data_d;     // p_pos

                const OPP_INT* opp_k1_dat0_stride_sycl = opp_k1_dat0_stride_s;
                const OPP_INT* opp_k1_dat1_stride_sycl = opp_k1_dat1_stride_s;

                const OPP_REAL* CONST_extents_sycl = CONST_extents_s;
                const OPP_REAL* CONST_dt_sycl = CONST_dt_s;

                // user provided elemental kernel
                // -----------------------------------------------------------------------------------------
                auto update_pos_kernel_sycl = [=](const double* part_vel, double* part_pos)
                {
                    for (int dm = 0; dm < 2; dm++) {

                        part_pos[(dm) * opp_k1_dat1_stride_sycl[0]] += 
                            part_vel[(dm) * opp_k1_dat0_stride_sycl[0]] * CONST_dt_sycl[0]; // s1 = s0 + ut

                        // correct for periodic boundary conditions
                        const int n_extent_offset_int =
                                sycl::fabs(part_pos[(dm)*opp_k1_dat1_stride_sycl[0]]) + 2.0;
                        const double temp_pos = part_pos[(dm) * opp_k1_dat1_stride_sycl[0]] + 
                                                n_extent_offset_int * CONST_extents_sycl[dm];
                        part_pos[(dm)*opp_k1_dat1_stride_sycl[0]] =
                                sycl::fmod((double)temp_pos, CONST_extents_sycl[dm]);
                    }
                };

                // -----------------------------------------------------------------------------------------
                auto kernel = [=](sycl::nd_item<1> item) {
                    const int tid = item.get_global_linear_id();
                    const int n = tid + start;
                    if (n < end) {
                        update_pos_kernel_sycl(
                            &arg0_data_d[n],
                            &arg1_data_d[n]);
                    }
                };

                // -----------------------------------------------------------------------------------------
                cgh.parallel_for(sycl::nd_range<1>(block_size*num_blocks,block_size), kernel);
            });
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();

    opp_profiler->end("update_pos_kernel");
}
