#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;

OPP_INT* opp_k3_dat0_stride_s = nullptr;
OPP_INT* opp_k3_dat1_stride_s = nullptr;

enum Dim {
    x = 0,
    y = 1,
};

//--------------------------------------------------------------
void opp_par_loop_all__verify_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // c_idx | OPP_READ
    opp_arg arg2 // | OPP_INC
)
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("verify_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__verify_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    OPP_INT *arg2_host_data = (OPP_INT *)args[2].data;

    opp_set_stride(opp_k3_dat0_stride_s, opp_k3_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat1_stride_s, opp_k3_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;

    int max_blocks = num_blocks;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_INT));
    reduction_size = MAX(reduction_size, sizeof(OPP_INT));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[2].data   = OPP_reduct_h + reduction_bytes;
    args[2].data_d = OPP_reduct_d + reduction_bytes;
    OPP_INT arg2_offset = reduction_bytes / sizeof(OPP_INT);

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_INT *)args[2].data)[b * 1 + d] = OPP_INT_ZERO;
    }

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_INT));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) {
        
        opp_queue->submit([&](sycl::handler &cgh) {

            auto opp_k3_dat0_stride_sycl = opp_k3_dat0_stride_s;
            auto CONST_cell_width_sycl = CONST_cell_width_s;
            auto CONST_ndimcells_sycl = CONST_ndimcells_s;

            const OPP_REAL *arg0_data_sycl = (OPP_REAL *)args[0].data_d;
            const OPP_INT *arg1_data_sycl = (OPP_INT *)args[1].data_d;
            const OPP_INT *p2c_map_sycl = (OPP_INT *)set->mesh_relation_dat->data_d;
            OPP_INT *arg2_data_sycl = (OPP_INT *)args[2].data_d;

            const OPP_INT set_size = iter_size;

            sycl::accessor<OPP_INT, 1, sycl::access::mode::read_write,
                    sycl::access::target::local> red_int(block_size, cgh); // temp

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto verify_kernel = [=](
                    const double* p_pos,
                    const int* c_gbl_idx,
                    int* incorrect_count) {

                // get the cell boundaries for the current cell_index - using global cell index
                int ix = -1, iy = -1;
                int _ix, _iy; _ix  = ((*c_gbl_idx)); _iy  = _ix/int(CONST_ndimcells_sycl[Dim::x]); _ix -= _iy*int(CONST_ndimcells_sycl[Dim::x]); (ix) = _ix; (iy) = _iy;;

                if (ix < 0 || iy < 0) {
                    // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
                    //     ix, iy, (*c_gbl_idx), CONST_ndimcells[Dim::x]);
                    (*incorrect_count)++;
                    return;
                }

                // get the boundaries of that cell
                const double boundary_ll[2] = { (ix * CONST_cell_width_sycl[0]), (iy * CONST_cell_width_sycl[0]) };

                // check whether the current particle is within those boundaries or not!
                const double p_pos_x = p_pos[(Dim::x) * opp_k3_dat0_stride_sycl[0]];
                if (p_pos_x < boundary_ll[Dim::x] ||
                        p_pos_x > (boundary_ll[Dim::x] + CONST_cell_width_sycl[0])) {

                    (*incorrect_count)++;
                    return;
                }

                const double p_pos_y = p_pos[(Dim::y) * opp_k3_dat0_stride_sycl[0]];
                if (p_pos_y < boundary_ll[Dim::y] ||
                        p_pos_y > (boundary_ll[Dim::y] + CONST_cell_width_sycl[0])) {

                    (*incorrect_count)++;
                    return;
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                OPP_INT gbl2_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl2_local[d] = OPP_INT_ZERO;

                for (int n = item.get_global_linear_id(); n < set_size; 
                                        n += item.get_global_range()[0]) {

                    const OPP_INT* opp_p2c = p2c_map_sycl + n;

                    verify_kernel(
                        arg0_data_sycl + n,
                        arg1_data_sycl + opp_p2c[0],
                        gbl2_local);
                }

                for (int d = 0; d < 1; d++){
                    opp_reduction<OPP_INC,1>(arg2_data_sycl, (arg2_offset + d + item.get_group_linear_id() * 1), 
                                    gbl2_local[d], red_int, item);
                }                    
            };

            // -----------------------------------------------------------------------------------------
            cgh.parallel_for(sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg2_host_data[d] += ((OPP_INT *)args[2].data)[b * 1 + d];
    }

    args[2].data = (char *)arg2_host_data;
    opp_mpi_reduce(&args[2], arg2_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();

    opp_profiler->end("verify_kernel");
}
