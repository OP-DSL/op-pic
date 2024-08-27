#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>

//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_c2c_map_stride = -1;

static dpct::constant_memory<OPP_INT, 0> opp_k2_dat0_stride_d;
static dpct::constant_memory<OPP_INT, 0> opp_k2_dat1_stride_d;
static dpct::constant_memory<OPP_INT, 0> opp_k2_c2c_map_stride_d;

namespace opp_k2 {
enum Dim {
    x = 0,
    y = 1,
};

enum CellMap {
    xd_y = 0,
    xu_y,
    x_yd,
    x_yu
};

inline void move_kernel(char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double* part_pos, const double* cell_pos_ll,
    OPP_INT opp_k2_dat0_stride_d, OPP_INT opp_k2_dat1_stride_d,
    OPP_INT opp_k2_c2c_map_stride_d, OPP_REAL const *CONST_cell_width_d)
{
    // check for x direction movement
    const double part_pos_x = part_pos[(Dim::x) * opp_k2_dat0_stride_d];
    if (part_pos_x < cell_pos_ll[(Dim::x) * opp_k2_dat1_stride_d]) {
        opp_p2c[0] = opp_c2c[(CellMap::xd_y) * opp_k2_c2c_map_stride_d];

        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (part_pos_x > (cell_pos_ll[(Dim::x) * opp_k2_dat1_stride_d] + CONST_cell_width_d[0])) {
        opp_p2c[0] = opp_c2c[(CellMap::xu_y) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    // check for y direction movement
    const double part_pos_y = part_pos[(Dim::y) * opp_k2_dat0_stride_d];
    if (part_pos_y < cell_pos_ll[(Dim::y) * opp_k2_dat1_stride_d]) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yd) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (part_pos_y > (cell_pos_ll[(Dim::y) * opp_k2_dat1_stride_d] + CONST_cell_width_d[0])) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yu) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    { opp_move_status_flag = OPP_MOVE_DONE; };
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
inline bool opp_part_check_status_cuda(char& move_flag, bool& iter_one_flag, 
        int* cell_id, int particle_index, int& remove_count, int *remove_particle_indices, 
        int *move_particle_indices, int *move_cell_indices, int *move_count,
        int OPP_cells_set_size_d) 
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE)
    {
        *cell_id = MAX_CELL_INDEX;
        const int removeArrayIndex =
            dpct::atomic_fetch_add<sycl::access::address_space::generic_space>(
                &remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (*cell_id >= OPP_cells_set_size_d)
    {
        // cell_id cell is not owned by the current mpi rank, need to communicate
        int moveArrayIndex =
            dpct::atomic_fetch_add<sycl::access::address_space::generic_space>(
                move_count, 1);
        move_particle_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *cell_id;

        // Needs to be removed from the current rank, 
        // particle packing will be done just prior exchange and removal
        move_flag = OPP_NEED_REMOVE;
        const int removeArrayIndex =
            dpct::atomic_fetch_add<sycl::access::address_space::generic_space>(
                &remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // cell_id is an own cell and move_flag == OPP_NEED_MOVE
    return true;
}
}

void opp_dev_move_kernel(
    const OPP_REAL *__restrict__ dat0,     // p_pos
    const OPP_REAL *__restrict__ dat1,     // c_pos_ll
    OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ c2c_map,
    OPP_INT *__restrict__ particle_remove_count,
    OPP_INT *__restrict__ particle_remove_indices,
    OPP_INT *__restrict__ move_particle_indices,
    OPP_INT *__restrict__ move_cell_indices,
    OPP_INT *__restrict__ move_count,
    const OPP_INT start,
    const OPP_INT end
,
    const sycl::nd_item<3> &item_ct1,
    OPP_INT opp_k2_dat0_stride_d,
    OPP_INT opp_k2_dat1_stride_d,
    OPP_INT opp_k2_c2c_map_stride_d,
    OPP_REAL const *CONST_cell_width_d,
    int OPP_cells_set_size_d,
    int OPP_comm_iteration_d) 
{
    const int thread_id = item_ct1.get_local_id(2) +
                          item_ct1.get_group(2) * item_ct1.get_local_range(2);

    if (thread_id + start < end) {

        const int n = thread_id + start;

        OPP_INT *opp_p2c = (p2c_map + n);
        char move_flag = OPP_NEED_MOVE;
        bool iter_one_flag = (OPP_comm_iteration_d > 0) ? false : true;

        do
        {
            const OPP_INT p2c = opp_p2c[0]; // get the value here, since the kernel might change it
            const OPP_INT* opp_c2c = c2c_map + p2c;

            opp_k2::move_kernel(move_flag, iter_one_flag, opp_c2c, opp_p2c,
                                dat0 + n, // p_pos
                                dat1 + p2c, opp_k2_dat0_stride_d,
                                opp_k2_dat1_stride_d, opp_k2_c2c_map_stride_d,
                                CONST_cell_width_d // c_pos_ll

            );

        } while (opp_k2::opp_part_check_status_cuda(
            move_flag, iter_one_flag, opp_p2c, n, *particle_remove_count,
            particle_remove_indices, move_particle_indices, move_cell_indices,
            move_count, OPP_cells_set_size_d));
    }
}

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_pos | OPP_READ
    opp_arg arg1   // c_pos_ll | OPP_READ
)
{
    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.in_order_queue();
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    if (OPP_cells_set_size != set->cells_set->size) {
        OPP_cells_set_size = set->cells_set->size;
        cutilSafeCall(
            DPCT_CHECK_ERROR(q_ct1
                                 .memcpy(OPP_cells_set_size_d.get_ptr(),
                                         &OPP_cells_set_size, sizeof(int))
                                 .wait()));
    }
    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;
    if (opp_k2_c2c_map_stride != c2c_stride) {
        opp_k2_c2c_map_stride = c2c_stride;
        cutilSafeCall(DPCT_CHECK_ERROR(
            q_ct1
                .memcpy(opp_k2_c2c_map_stride_d.get_ptr(),
                        &opp_k2_c2c_map_stride, sizeof(OPP_INT))
                .wait()));
    }

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    int num_blocks = 200;

    do 
    {
        if (opp_k2_dat0_stride != args[0].dat->set->set_capacity) {
            opp_k2_dat0_stride = args[0].dat->set->set_capacity;
            cutilSafeCall(DPCT_CHECK_ERROR(
                q_ct1
                    .memcpy(opp_k2_dat0_stride_d.get_ptr(), &opp_k2_dat0_stride,
                            sizeof(OPP_INT))
                    .wait()));
        }
        if (opp_k2_dat1_stride != args[1].dat->set->set_capacity) {
            opp_k2_dat1_stride = args[1].dat->set->set_capacity;
            cutilSafeCall(DPCT_CHECK_ERROR(
                q_ct1
                    .memcpy(opp_k2_dat1_stride_d.get_ptr(), &opp_k2_dat1_stride,
                            sizeof(OPP_INT))
                    .wait()));
        }

        opp_init_particle_move(set, nargs, args);
        cutilSafeCall(
            DPCT_CHECK_ERROR(q_ct1
                                 .memcpy(OPP_comm_iteration_d.get_ptr(),
                                         &OPP_comm_iteration, sizeof(int))
                                 .wait()));

        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

        /*
        DPCT1049:1: The work-group size passed to the SYCL kernel may exceed the
        limit. To get the device limit, query info::device::max_work_group_size.
        Adjust the work-group size if needed.
        */
        {
            extern dpct::constant_memory<OPP_REAL, 1> CONST_cell_width_d;
            extern dpct::constant_memory<int, 0> OPP_cells_set_size_d;
            extern dpct::constant_memory<int, 0> OPP_comm_iteration_d;

            opp_k2_dat0_stride_d.init();
            opp_k2_dat1_stride_d.init();
            opp_k2_c2c_map_stride_d.init();
            CONST_cell_width_d.init();
            OPP_cells_set_size_d.init();
            OPP_comm_iteration_d.init();

            q_ct1.submit([&](sycl::handler &cgh) {
                auto opp_k2_dat0_stride_d_ptr_ct1 =
                    opp_k2_dat0_stride_d.get_ptr();
                auto opp_k2_dat1_stride_d_ptr_ct1 =
                    opp_k2_dat1_stride_d.get_ptr();
                auto opp_k2_c2c_map_stride_d_ptr_ct1 =
                    opp_k2_c2c_map_stride_d.get_ptr();
                auto CONST_cell_width_d_ptr_ct1 = CONST_cell_width_d.get_ptr();
                auto OPP_cells_set_size_d_ptr_ct1 =
                    OPP_cells_set_size_d.get_ptr();
                auto OPP_comm_iteration_d_ptr_ct1 =
                    OPP_comm_iteration_d.get_ptr();

                const double *args_data_d_ct0 = (OPP_REAL *)args[0].data_d;
                const double *args_data_d_ct1 = (OPP_REAL *)args[1].data_d;
                int *args_data_d_ct2 = (OPP_INT *)args[2].data_d;
                const int *c2c_map_map_d_ct3 = (OPP_INT *)c2c_map->map_d;
                int *set_particle_remove_count_d_ct4 =
                    (OPP_INT *)set->particle_remove_count_d;
                int *OPP_remove_particle_indices_d_ct5 =
                    (OPP_INT *)OPP_remove_particle_indices_d;
                int *OPP_move_particle_indices_d_ct6 =
                    (OPP_INT *)OPP_move_particle_indices_d;
                int *OPP_move_cell_indices_d_ct7 =
                    (OPP_INT *)OPP_move_cell_indices_d;
                int *OPP_move_count_d_ct8 = (OPP_INT *)OPP_move_count_d;
                const int OPP_iter_start_ct9 = OPP_iter_start;
                const int OPP_iter_end_ct10 = OPP_iter_end;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, num_blocks) *
                                          sycl::range<3>(1, 1, block_size),
                                      sycl::range<3>(1, 1, block_size)),
                    [=](sycl::nd_item<3> item_ct1) {
                        opp_dev_move_kernel(
                            args_data_d_ct0, args_data_d_ct1, args_data_d_ct2,
                            c2c_map_map_d_ct3, set_particle_remove_count_d_ct4,
                            OPP_remove_particle_indices_d_ct5,
                            OPP_move_particle_indices_d_ct6,
                            OPP_move_cell_indices_d_ct7, OPP_move_count_d_ct8,
                            OPP_iter_start_ct9, OPP_iter_end_ct10, item_ct1,
                            *opp_k2_dat0_stride_d_ptr_ct1,
                            *opp_k2_dat1_stride_d_ptr_ct1,
                            *opp_k2_c2c_map_stride_d_ptr_ct1,
                            CONST_cell_width_d_ptr_ct1,
                            *OPP_cells_set_size_d_ptr_ct1,
                            *OPP_comm_iteration_d_ptr_ct1);
                    });
            });
        }

    } while (opp_finalize_particle_move(set)); 

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(DPCT_CHECK_ERROR(dev_ct1.queues_wait_and_throw()));

    opp_profiler->end("move_kernel");
}
void opp_init_direct_hop_cg(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1 // c_pos_ll | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    opp_profiler->end("Setup_Mover");
}
