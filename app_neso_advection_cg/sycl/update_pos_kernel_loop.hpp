#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;

static dpct::constant_memory<OPP_INT, 0> opp_k1_dat0_stride_d;
static dpct::constant_memory<OPP_INT, 0> opp_k1_dat1_stride_d;

namespace opp_k1 {
inline void update_pos_kernel(const double* part_vel, double* part_pos,
                              OPP_INT opp_k1_dat0_stride_d,
                              OPP_INT opp_k1_dat1_stride_d,
                              OPP_REAL const *CONST_extents_d,
                              OPP_REAL const *CONST_dt_d)
{
    for (int dm = 0; dm < 2; dm++) {

        part_pos[(dm) * opp_k1_dat1_stride_d] += part_vel[(dm) * opp_k1_dat0_stride_d] * CONST_dt_d[0]; // s1 = s0 + ut

        // correct for periodic boundary conditions
        const int n_extent_offset_int =
            sycl::fabs(part_pos[(dm)*opp_k1_dat1_stride_d]) + 2.0;
        const double temp_pos = part_pos[(dm) * opp_k1_dat1_stride_d] + n_extent_offset_int * CONST_extents_d[dm];
        part_pos[(dm)*opp_k1_dat1_stride_d] =
            sycl::fmod((double)temp_pos, CONST_extents_d[dm]);
    }
}

}

//--------------------------------------------------------------
void opp_dev_update_pos_kernel(
    const OPP_REAL *__restrict__ dat0,  // p_vel
    OPP_REAL *__restrict__ dat1,  // p_pos
    const OPP_INT start,
    const OPP_INT end
,
    const sycl::nd_item<3> &item_ct1,
    OPP_INT opp_k1_dat0_stride_d,
    OPP_INT opp_k1_dat1_stride_d,
    OPP_REAL const *CONST_extents_d,
    OPP_REAL const *CONST_dt_d) 
{
    const int thread_id = item_ct1.get_local_id(2) +
                          item_ct1.get_group(2) * item_ct1.get_local_range(2);

    for (int n = thread_id; n < (end - start);
         n += item_ct1.get_local_range(2) * item_ct1.get_group_range(2)) {

        opp_k1::update_pos_kernel(dat0 + n, // p_vel
                                  dat1 + n, opp_k1_dat0_stride_d,
                                  opp_k1_dat1_stride_d, CONST_extents_d,
                                  CONST_dt_d // p_pos
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__update_pos_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1 // p_pos | OPP_RW
)
{
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.in_order_queue();
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    if (opp_k1_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k1_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(
            DPCT_CHECK_ERROR(q_ct1
                                 .memcpy(opp_k1_dat0_stride_d.get_ptr(),
                                         &opp_k1_dat0_stride, sizeof(OPP_INT))
                                 .wait()));
    }
    if (opp_k1_dat1_stride != args[1].dat->set->set_capacity) {
        opp_k1_dat1_stride = args[1].dat->set->set_capacity;
        cutilSafeCall(
            DPCT_CHECK_ERROR(q_ct1
                                 .memcpy(opp_k1_dat1_stride_d.get_ptr(),
                                         &opp_k1_dat1_stride, sizeof(OPP_INT))
                                 .wait()));
    }

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
            /*
            DPCT1049:0: The work-group size passed to the SYCL kernel may exceed
            the limit. To get the device limit, query
            info::device::max_work_group_size. Adjust the work-group size if
            needed.
            */
            extern dpct::constant_memory<OPP_REAL, 1> CONST_extents_d;
            extern dpct::constant_memory<OPP_REAL, 1> CONST_dt_d;

            opp_k1_dat0_stride_d.init();
            opp_k1_dat1_stride_d.init();
            CONST_extents_d.init();
            CONST_dt_d.init();

            dpct::has_capability_or_fail(q_ct1.get_device(),
                                         {sycl::aspect::fp64});

            q_ct1.submit([&](sycl::handler &cgh) {
                auto opp_k1_dat0_stride_d_ptr_ct1 =
                    opp_k1_dat0_stride_d.get_ptr();
                auto opp_k1_dat1_stride_d_ptr_ct1 =
                    opp_k1_dat1_stride_d.get_ptr();
                auto CONST_extents_d_ptr_ct1 = CONST_extents_d.get_ptr();
                auto CONST_dt_d_ptr_ct1 = CONST_dt_d.get_ptr();

                const double *args_data_d_ct0 = (OPP_REAL *)args[0].data_d;
                double *args_data_d_ct1 = (OPP_REAL *)args[1].data_d;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, num_blocks) *
                                          sycl::range<3>(1, 1, block_size),
                                      sycl::range<3>(1, 1, block_size)),
                    [=](sycl::nd_item<3> item_ct1) {
                        opp_dev_update_pos_kernel(
                            args_data_d_ct0, args_data_d_ct1, start, end,
                            item_ct1, *opp_k1_dat0_stride_d_ptr_ct1,
                            *opp_k1_dat1_stride_d_ptr_ct1,
                            CONST_extents_d_ptr_ct1, CONST_dt_d_ptr_ct1);
                    });
            });
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(DPCT_CHECK_ERROR(dev_ct1.queues_wait_and_throw()));

    opp_profiler->end("update_pos_kernel");
}
