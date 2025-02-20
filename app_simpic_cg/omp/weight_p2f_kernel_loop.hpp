
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k3 {
void weight_p2f_kernel(
        double* node0_field_J,
        double* node1_field_J,
        const double* particle0_position_x
    )
{
    double xx = ((particle0_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    (*node0_field_J) += (CONST_qscale[0] * (1.0 - frac));  // Can change qscale to be calculated from particle data
    (*node1_field_J) += (CONST_qscale[0] * frac);
}
}

void opp_par_loop_all__weight_p2f_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_j | OPP_INC
    opp_arg arg1, // n_field_j | OPP_INC
    opp_arg arg2 // p_pos_x | OPP_READ
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("weight_p2f_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_p2f_kernel set_size %d", set->size);

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
#ifdef USE_MPI 
    opp_init_double_indirect_reductions(nargs, args);
#endif
 
    opp_mpi_halo_wait_all(nargs, args);    
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 

    opp_create_thread_level_data<OPP_REAL>(args[0], 0);
  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;

        OPP_REAL* arg0_dat_thread_data = (OPP_REAL *)((*(args[0].dat->thread_data))[thr]);
      
        for (size_t n = start; n < finish; n++)
        { 
            OPP_INT* opp_p2c = OPP_mesh_relation_data + n;
            const OPP_INT *map0 = args[0].map_data + (opp_p2c[0] * 2);
   
            opp_k3::weight_p2f_kernel(
                arg0_dat_thread_data + (map0[0] * 1),
                arg0_dat_thread_data + (map0[1] * 1),
                (const OPP_REAL *)args[2].data + (n * 1)
            );
        }
    }
    
    opp_reduce_thread_level_data<OPP_REAL>(args[0]);

#ifdef USE_MPI     
    opp_exchange_double_indirect_reductions(nargs, args);
    opp_complete_double_indirect_reductions(nargs, args);
#endif
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("weight_p2f_kernel");
}
