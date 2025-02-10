
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

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    

    OPP_REAL arg1_l[nthreads * 1];
    for (int thr = 0; thr < nthreads; thr++)
        for (int d = 0; d < 1; d++)
            arg1_l[1 * thr + d] = OPP_REAL_ZERO;

  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            opp_k8::get_sigma_ef_sq_kernel(
                (const OPP_REAL *)args[0].data + (n * 3),
                arg1_l + (1 * thr)
            );
        }
    }
    
    for (int thr = 0; thr < nthreads; thr++) 
        for (int d = 0; d < 1; d++)
            ((OPP_REAL*)args[1].data)[d] += (arg1_l[1 * thr + d]);
#ifdef USE_MPI 
    opp_mpi_reduce(&args[1], (OPP_REAL *)args[1].data);
#endif
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("get_sigma_ef_sq_kernel");
}
