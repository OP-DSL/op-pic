

inline void calculate_injection_distribution(
    int* injected_total,
    double* face_area,
    int* particle_distribution,
    double* remainder 
)
{   
    /*number of real ions per sec, given prescribed density and velocity*/
    double num_per_sec = CONST_plasma_den * CONST_ion_velocity * (*face_area);

    /*number of ions to generate in this time step*/
    double num_real = num_per_sec * CONST_dt;

    /*fraction number of macroparticles*/
    double fnum_mp = num_real / CONST_spwt + (*remainder);

    /*integer number of macroparticles*/
    int num_mp = (int)fnum_mp;

    /*update reminder*/
    (*remainder) = fnum_mp - num_mp;

    (*injected_total) += num_mp;

    (*particle_distribution) = (*injected_total);
}

void oppic_seq_loop_inject__Increase_particle_count(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
)
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_seq_loop_inject__Increase_particle_count num_particles %d diff %d\n", set->size, set->diff);

    int nargs = 5; // Add one more to have mesh_relation arg
    oppic_arg args[nargs] = { arg0, arg1, arg2, arg3, oppic_arg_dat(particles_set->cell_index_dat, OP_RW) };
    particles_set->cell_index_dat->dirty_hd = Dirty::Host; // make mesh relation dirty and download new data from device

    int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, Device_CPU);

    for (int i = 0; i < set_size; i++)
    {   
        calculate_injection_distribution(
            ((int *)arg0.data),
            &((double *)arg1.data)[i],
            &((int *)arg2.data)[i],
            ((double *)arg3.data) 
        );
    }

    oppic_increase_particle_count(particles_set, *((int *)arg0.data));

    int* part_mesh_connectivity = (int *)particles_set->cell_index_dat->data;
    int* distribution           = (int *)arg2.data;

    int start = (particles_set->size - particles_set->diff);
    int j = 0;

    for (int i = 0; i < particles_set->diff; i++)
    {
        if (i >= distribution[j]) j++; // check whether it is j or j-1    
        part_mesh_connectivity[start + i] = j;
    }  

    op_mpi_set_dirtybit(nargs, args);
}


