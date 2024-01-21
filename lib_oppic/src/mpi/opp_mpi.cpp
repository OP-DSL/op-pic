
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

#include <opp_mpi.h>

opp_move_var move_var;
void opp_part_pack(oppic_set set);
void opp_part_unpack(oppic_set set);


//*******************************************************************************
void opp_init(int argc, char **argv) 
{

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULL, "opp::PetscMPI");
#else
    MPI_Init(&argc, &argv);
#endif

    OP_MPI_WORLD = MPI_COMM_WORLD;
    OP_MPI_GLOBAL = MPI_COMM_WORLD;
    
    MPI_Comm_rank(OP_MPI_WORLD, &OPP_rank);
    MPI_Comm_size(OP_MPI_WORLD, &OPP_comm_size);

    if (OPP_rank == OPP_ROOT) 
        opp_printf("oppic_init", "initialized with %d MPI ranks", OPP_comm_size);
    
    oppic_init_core(argc, argv);

    opp_profiler->reg("Mv_Finalize");
}

//*******************************************************************************
void opp_exit() 
{
    if (OP_DEBUG) opp_printf("oppic_exit", "");

    globalMover.reset();
    cellMapper.reset();
    boundingBox.reset();
    comm.reset();
    
    {   
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
        
        // free(set_import_buffer_size);

        // for (int i = 0; i < OP_import_index; i++)
        //     free(OP_import_list[i]);
        // if (OP_import_list)
        //     free(OP_import_list);
        
        // for (int i = 0; i < OP_export_index; i++)
        //     free(OP_export_list[i]);
        // if (OP_export_list)
        //     free(OP_export_list);
    }

    oppic_exit_core();

#ifdef USE_PETSC
    PetscFinalize();
#else
    MPI_Finalize();
#endif
    
}

//****************************************
void opp_abort(std::string s)
{
    opp_printf("opp_abort", "%s", s.c_str());
    MPI_Abort(OP_MPI_WORLD, 2);
}

//****************************************
oppic_set opp_decl_mesh_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map opp_decl_mesh_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    oppic_map map = oppic_decl_map_core(from, to, dim, imap, name);

    if (OP_DEBUG) opp_printf("oppic_decl_map", OPP_rank, " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

    return map;
}

//****************************************
oppic_dat opp_decl_mesh_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);
}

//****************************************
oppic_map oppic_decl_map_txt(oppic_set from, oppic_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)oppic_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    oppic_map map = opp_decl_mesh_map(from, to, dim, map_data, name);

    free(map_data);

    return map;
}

//****************************************
oppic_dat oppic_decl_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    oppic_dat dat = oppic_decl_dat_core(set, dim, type.c_str(), size, dat_data, name);

    free(dat_data);

    return dat;
}

//****************************************
oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, mapping);
}

//****************************************
oppic_arg opp_get_arg(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, acc, mapping);
}
oppic_arg opp_get_arg(oppic_dat dat, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, acc, mapping);
}
oppic_arg opp_get_arg(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, acc, mapping);
}
oppic_arg opp_get_arg(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, idx, map, acc, mapping);
}


//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg opp_get_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg opp_get_arg_gbl(int *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg opp_get_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
oppic_set opp_decl_part_set(int size, char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(size, name, cells_set);
}
oppic_set opp_decl_part_set(char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(name, cells_set);
}

//****************************************
oppic_dat opp_decl_part_dat(oppic_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return oppic_decl_particle_dat_core(set, dim, type.c_str(), size, (char*)data, name, cell_index);
}

//****************************************
oppic_dat oppic_decl_particle_dat_txt(oppic_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)oppic_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    oppic_dat dat = oppic_decl_particle_dat_core(set, dim, type.c_str(), size, dat_data, name, cell_index);

    free(dat_data);

    return dat;
}

//****************************************
void oppic_increase_particle_count(oppic_set particles_set, const int num_particles_to_insert)
{
    if (!oppic_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_printf("oppic_increase_particle_count", "Error: oppic_increase_particle_count_core failed for particle set [%s]", particles_set->name);
        opp_abort("oppic_increase_particle_count");        
    }
}

//****************************************
void opp_inc_part_count_with_distribution(oppic_set particles_set, int num_particles_to_insert, oppic_dat part_dist, bool calc_new)
{
    if (OP_DEBUG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        opp_abort("opp_inc_part_count_with_distribution_core");        
    }
}

//****************************************
void oppic_reset_num_particles_to_insert(oppic_set set) // unused
{
    oppic_reset_num_particles_to_insert_core(set);
}

//****************************************
void oppic_particle_sort(oppic_set set) // unused
{ 
    oppic_particle_sort_core(set);
}

//****************************************
void opp_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_m" + std::to_string(OPP_comm_size);
    oppic_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    oppic_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

//****************************************
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args)
{ 

    oppic_init_particle_move_core(set);

    move_part_indices.clear();

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }
    else
    {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++)
        {
            if (args[i].argtype == OP_ARG_DAT && args[i].dat->set->is_particle)
            {
                args[i].data = args[i].dat->data;
            }
        }
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
}

//****************************************
bool opp_finalize_particle_move(oppic_set set)
{ 

    if (OP_DEBUG) opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

    opp_process_marked_particles(set); 

    opp_part_pack(set);
    
    // send the counts and send the particles  
    opp_part_exchange(set);  

    // Can fill the holes here, since the communicated particles will be added at the end
    oppic_finalize_particle_move_core(set);

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) opp_printf("opp_finalize_particle_move", "auto sorting particle set [%s]", set->name);
        oppic_particle_sort(set);
    }

    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");

        return false; // all mpi ranks do not have anything to communicate to any rank
    }
        
    opp_part_wait_all(set); // wait till all the particles are communicated
    
    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_part_unpack(set);

    // cleanSendRecvBuffers(set);

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");

    return true;
}

//****************************************
void opp_reset_dat(oppic_dat dat, char* val, opp_reset reset)
{
    if (!val)
    {
        opp_printf("opp_reset_dat", "Error: val is NULL");
        return;
    }

    int start = -1;
    int end = -1;

    opp_get_start_end(dat->set, reset, start, end);

    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }

    // TODO : Check whether this is OK for all the reset options!
    dat->dirtybit = 0;
}

//*******************************************************************************
void opp_partition(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    opp_profiler->start("opp_partition");

    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();

    MPI_Barrier(MPI_COMM_WORLD);
    
    opp_profiler->end("opp_partition");
}

//*******************************************************************************
void opp_mpi_print_dat_to_txtfile(op_dat dat, const char *file_name) 
{
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + std::to_string(OPP_comm_size) + std::string("_") + file_name;
    // rearrange data back to original order in mpi
    opp_dat temp = opp_mpi_get_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    free(temp->data);
    free(temp->set);
    free(temp);
}

//*******************************************************************************
opp_move_var opp_get_move_var(int thread)
{
    // no perf improvement by using a buffered move var, could create a new here instead
    
    move_var.move_status = OPP_MOVE_DONE;

    if (OPP_comm_iteration != 0) // TRUE means communicated particles, no need to do the iteration one calculations
        move_var.iteration_one = false;
    else
        move_var.iteration_one = true;

    return move_var; // passing the object for now :(
}

//*******************************************************************************
opp_dat opp_fetch_data(opp_dat dat) {
    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    // rearrange data backe to original order in mpi
    return opp_mpi_get_data(dat);
}