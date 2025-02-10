
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k8 {
inline void get_sigma_ef_sq_kernel(
    const double* val,
    double* max_val
) {
    for (int dim = 0; dim < 3; ++dim) {
        *max_val += (val[dim] * val[dim]);
    }
}
}

void opp_par_loop_all__get_sigma_ef_sq_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1 // | OPP_INC
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("get_sigma_ef_sq_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_sigma_ef_sq_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    for (int n = 0; n < iter_size; ++n) 
    {
        if (n == set->core_size)
            opp_mpi_halo_wait_all(nargs, args);

        opp_k8::get_sigma_ef_sq_kernel(
            (const OPP_REAL *)args[0].data + (n * 3),
            (OPP_REAL *)args[1].data
        );
    }
    if (iter_size == 0 || iter_size == set->core_size)
        opp_mpi_halo_wait_all(nargs, args);
    opp_mpi_reduce(&args[1], (OPP_REAL *)args[1].data);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("get_sigma_ef_sq_kernel");
}
