
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_dat2_stride = -1;

__constant__ OPP_INT opp_k3_dat0_stride_d;
__constant__ OPP_INT opp_k3_dat1_stride_d;
__constant__ OPP_INT opp_k3_dat2_stride_d;



namespace opp_k3 {
__device__ inline void calculate_new_pos_vel_kernel(
    const double *cell_ef,
    double *part_pos,
    double *part_vel ) {

    const double coefficient1 = CONST_charge_d[0] / CONST_mass_d[0] * (CONST_dt_d[0]);
    for (int i = 0; i < 3; i++) {
        part_vel[(i) * opp_k3_dat2_stride_d] += (coefficient1 * cell_ef[(i) * opp_k3_dat0_stride_d]);
        part_pos[(i) * opp_k3_dat1_stride_d] += part_vel[(i) * opp_k3_dat2_stride_d] * (CONST_dt_d[0]);
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_calculate_new_pos_vel_kernel(
    const OPP_REAL *__restrict__ dat0,  // c_ef
    OPP_REAL *__restrict__ dat1,  // p_pos
    OPP_REAL *__restrict__ dat2,  // p_vel
    const OPP_INT *__restrict__ p2c_map,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        const OPP_INT* opp_p2c = p2c_map + n;
        
        opp_k3::calculate_new_pos_vel_kernel(
            dat0 + opp_p2c[0], // c_ef 
            dat1 + n, // p_pos 
            dat2 + n // p_vel 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__calculate_new_pos_vel_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1, // p_pos | OPP_WRITE
    opp_arg arg2 // p_vel | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("calculate_new_pos_vel_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__calculate_new_pos_vel_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    if (opp_k3_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k3_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k3_dat0_stride_d, &opp_k3_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k3_dat1_stride != args[1].dat->set->set_capacity) {
        opp_k3_dat1_stride = args[1].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k3_dat1_stride_d, &opp_k3_dat1_stride, sizeof(OPP_INT)));
    }
    if (opp_k3_dat2_stride != args[2].dat->set->set_capacity) {
        opp_k3_dat2_stride = args[2].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k3_dat2_stride_d, &opp_k3_dat2_stride, sizeof(OPP_INT)));
    }

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
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
            opp_dev_calculate_new_pos_vel_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_ef
                (OPP_REAL *)args[1].data_d,     // p_pos
                (OPP_REAL *)args[2].data_d,     // p_vel
                (OPP_INT *)set->mesh_relation_dat->data_d,
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());   
 
    opp_profiler->end("calculate_new_pos_vel_kernel");
}
