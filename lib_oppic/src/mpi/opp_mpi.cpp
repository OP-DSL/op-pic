
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


MPI_Comm OP_MPI_WORLD;
MPI_Comm OP_MPI_GLOBAL;

opp_move_var move_var;

op_dat op_mpi_get_data(op_dat dat);

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

    if (OP_DEBUG) opp_printf("oppic_init", "");
    
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
void oppic_reset_num_particles_to_insert(oppic_set set)
{
    oppic_reset_num_particles_to_insert_core(set);
}

// //****************************************
// void oppic_mark_particle_to_remove(oppic_set set, int particle_index)
// {
//     oppic_mark_particle_to_remove_core(set, particle_index);
// }

// //****************************************
// void oppic_remove_marked_particles_from_set(oppic_set set)
// {
//     oppic_remove_marked_particles_from_set_core(set);
// }
// void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove)
// {
//     oppic_remove_marked_particles_from_set_core(set, idx_to_remove);
// }

//****************************************
void oppic_particle_sort(oppic_set set)
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

// //****************************************
// void oppic_dump_dat(oppic_dat dat)
// {
//     oppic_dump_dat_core(dat);
// }

//****************************************
void opp_init_particle_move(oppic_set set, int nargs, oppic_arg *args)
{ 

    oppic_init_particle_move_core(set);

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

// //****************************************
// void oppic_mark_particle_to_move(oppic_set set, int particle_index, int move_status)
// {
//     oppic_mark_particle_to_move_core(set, particle_index, move_status);
// }

//****************************************
bool opp_finalize_particle_move(oppic_set set)
{ 

    if (OP_DEBUG) opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

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
        
    opp_part_wait_all(set); // wait till all the particles are communicated and added to the dats

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

    switch (reset)
    {
        case OPP_Reset_Core:
            start = 0;
            end = dat->set->core_size;
            break;
        case OPP_Reset_Set:
            start = 0;
            end = dat->set->size;
            break;
        case OPP_Reset_All:
            start = 0;
            end = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
            break;
        case OPP_Reset_ieh:
            start = dat->set->size;
            end = dat->set->size + dat->set->exec_size;
            break;
        case OPP_Reset_inh:
            start = dat->set->size + dat->set->exec_size;
            end = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
            break;
        default:
            opp_printf("opp_reset_dat", "Error: opp_reset failure");
    }

    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }

    // TODO : Check whether this is OK for all the reset options!
    dat->dirtybit = 0;
}

//****************************************
void opp_mpi_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        // TODO : Do not include double indirect reductions
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) && 
            (args[n].acc == OP_WRITE || args[n].acc == OP_RW || 
                (args[n].acc == OP_INC)))  //  && !is_double_indirect_reduction(args[n])
        {
            args[n].dat->dirtybit = 1;
            // args[n].dat->dirty_hd = Dirty::Device;
        }
    }
}

//*******************************************************************************
void opp_partition(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();
}

//*******************************************************************************
void opp_partition_core(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    if (lib_name == "PARMETIS_KWAY")
    {
#ifdef HAVE_PARMETIS
        if (prime_map != NULL)
        {
            opp_partition_kway(prime_map); // use parmetis kway partitioning
        }
        else
        {
            opp_abort("opp_partition PARMETIS_KWAY Error: Partitioning prime_map : NULL - UNSUPPORTED Partitioner Specification");  
        }
#else
        opp_abort("opp_partition_core PARMETIS_KWAY Error: Parmetis not installed or not defined");
#endif
    }
    else if (lib_name == "PARMETIS_GEOM")
    {
#ifdef HAVE_PARMETIS
        if (data != NULL)
        {
            opp_partition_geom(data); // use parmetis geometric partitioning
        }
        else
        {
            opp_abort("opp_partition PARMETIS_GEOM Error: Partitioning geom dat : NULL - UNSUPPORTED Partitioner Specification"); 
        }
#else
        opp_abort("opp_partition_core PARMETIS_GEOM Error: Parmetis not installed or not defined");
#endif
    }
    else if (lib_name == "EXTERNAL")
    {
        if (data != NULL)
        {
            opp_partition_external(prime_set, data); // use external partitioning dat
        }
        else
        {
            opp_abort("opp_partition EXTERNAL Error: Partitioning color dat : NULL - UNSUPPORTED Partitioner Specification"); 
        }
    }
    else if (lib_name != "")
    {
        opp_abort("opp_partition Error: Unsupported lib_name - UNSUPPORTED Partitioner Specification");
    }

    opp_halo_create();

    // TODO : I think this sanity check is wrong, in cell->nodes mapping there can be nodes indexed in import non-exec halo,
    // hence the mapping is obvously greater than to_set->size 
        // sanity check to identify if the partitioning results in ophan elements
        // int ctr = 0;
        // for (int i = 0; i < prime_map->from->size; i++)
        // {
        //     if (prime_map->map[2 * i] >= prime_map->to->size && prime_map->map[2 * i + 1] >= prime_map->to->size) 
        //         ctr++;
        // }
        // opp_printf("opp_partition()", "%s Orphans in prime map [%s]: %d", (ctr > 0) ? "Error:" : "", prime_map->name, ctr);

    opp_part_comm_init(); 

    std::vector<std::vector<int>> set_sizes(oppic_sets.size());

    for (oppic_set set : oppic_sets)
    {
        std::vector<int>& recv_vec = set_sizes[set->index];
        recv_vec.resize(OPP_comm_size * 3);

        std::vector<int> sizes{ set->size, set->exec_size, set->nonexec_size };
        MPI_Gather(&(sizes[0]), 3, MPI_INT, &(recv_vec[0]), 3, MPI_INT, OPP_ROOT, OP_MPI_WORLD);
    }

    // print the set sizes of all ranks after partitioning
    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "";

        for (oppic_set set : oppic_sets)
            log += "\t - " + std::string(set->name);

        opp_printf("opp_partition()", "(size|ieh|inh) %s", log.c_str());

        for (int i = 0; i < OPP_comm_size; i++)
        {
            log = "RANK [" + std::to_string(i) + "]";
            
            for (int j = 0; j < (int)oppic_sets.size(); j++)
                log += "\t- " + std::to_string(set_sizes[j][i * 3]) + "|" + 
                    std::to_string(set_sizes[j][i * 3 + 1]) + "|" + std::to_string(set_sizes[j][i * 3 + 2]);

            opp_printf("opp_partition()", "%s", log.c_str());
        }
    }
}

std::map<int, oppic_dat> negative_mapping_indices;

//*******************************************************************************
void opp_sanitize_all_maps()
{
    for (int i = 0; i < (int)oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_sanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        std::string name = std::string("AUTO_DAT_") + map->name;
        oppic_dat dat = opp_decl_mesh_dat(map->from, map->dim, DT_INT, (char*)map->map, name.c_str());  
        negative_mapping_indices[map->index] = dat;

        memset(dat->data, 0, (map->from->size * dat->size));

        for (int n = 0; n < map->from->size; n++)
        {
            int positive_mapping = -1;
            std::vector<int> index;

            for (int d = 0; d < map->dim; d++)
            {
                if (map->map[n * map->dim + d] < 0)
                {
                    index.push_back(n * map->dim + d);
                }
                else
                {
                    positive_mapping = map->map[n * map->dim + d];
                }
            }

            if (positive_mapping >= 0)
            {
                for (int i : index)
                {
                    map->map[i] = positive_mapping;
                    ((int*)dat->data)[i] = 1;
                }
            }
            else
            {
                opp_printf("opp_sanitize_all_maps", "Error: No positive mapping found at %d in map: %s", n, map->name);
            }
        }
    }

    // if (OP_DEBUG)
    // {
    //     for (int i = 0; i < oppic_maps.size(); i++)
    //     {
    //         oppic_map map = oppic_maps[i];
            
    //         opp_printf("opp_sanitize_all_maps", OPP_rank, " map: %s | from->size: %d | dim: %d", map->name, map->from->size, map->dim);

    //         for (int n = 0; n < map->from->size; n++)
    //         {
    //             for (int d = 1; d < map->dim; d++)
    //             {
    //                 if (map->map[n * map->dim + d] < 0)
    //                 {
    //                     opp_printf("opp_sanitize_all_maps", OPP_rank, "Error: map: %s | ptr: %p | negative mapping at index: %d [%d]", map->name, map->map, n, n * map->dim + d);
    //                 }
    //             }
    //         }
    //     }
    // }
}

//*******************************************************************************
void opp_desanitize_all_maps()
{
    for (int i = 0; i < (int)oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_desanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        auto it = negative_mapping_indices.find(map->index);
        if (it == negative_mapping_indices.end())
        {
            opp_printf("opp_desanitize_all_maps", "Error: Negative mappings not found for map: %s", map->name);
            continue;
        }
            
        oppic_dat dat = it->second;

        for (int x = 0; x < (map->from->size + map->from->exec_size) * map->dim; x++)
        {
            if (((int*)dat->data)[x] == 1)
                map->map[x] = -1;
        }
    }

    // could realloc the dat to a lower size, if required
    negative_mapping_indices.clear();
}


void opp_mpi_print_dat_to_txtfile(op_dat dat, const char *file_name) 
{
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + std::to_string(OPP_comm_size) + std::string("_") + file_name;
    // rearrange data back to original order in mpi
    op_dat temp = op_mpi_get_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    free(temp->data);
    free(temp->set);
    free(temp);
}

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
//*******************************************************************************

#define BOUNDING_TOLERENCE 1e-12

using namespace opp;

Comm::Comm(MPI_Comm comm_parent) {
    this->comm_parent = comm_parent;

    int rank_parent;
    CHECK(MPI_Comm_rank(comm_parent, &rank_parent))
    CHECK(MPI_Comm_split_type(comm_parent, MPI_COMM_TYPE_SHARED, 0,
                            MPI_INFO_NULL, &this->comm_intra))

    int rank_intra;
    CHECK(MPI_Comm_rank(this->comm_intra, &rank_intra))
    const int colour_intra = (rank_intra == 0) ? 1 : MPI_UNDEFINED;
    CHECK(MPI_Comm_split(comm_parent, colour_intra, 0, &this->comm_inter))

    CHECK(MPI_Comm_rank(this->comm_parent, &this->rank_parent))
    CHECK(MPI_Comm_rank(this->comm_intra, &this->rank_intra))
    CHECK(MPI_Comm_size(this->comm_parent, &this->size_parent))
    CHECK(MPI_Comm_size(this->comm_intra, &this->size_intra))
    
    if (comm_inter != MPI_COMM_NULL) {

        CHECK(MPI_Comm_rank(this->comm_inter, &this->rank_inter))
        CHECK(MPI_Comm_size(this->comm_inter, &this->size_inter))
    }

    if (OP_DEBUG)
        opp_printf("Comm", "rank_parent %d|s=%d rank_intra %d|s=%d rank_inter %d|s=%d",
            this->rank_parent, this->size_parent, this->rank_intra, this->size_intra, 
            this->rank_inter, this->size_inter);
};

Comm::~Comm() {

    if ((this->comm_intra != MPI_COMM_NULL) && (this->comm_intra != MPI_COMM_WORLD)) {
        
        CHECK(MPI_Comm_free(&this->comm_intra))
        this->comm_intra = MPI_COMM_NULL;
    }

    if ((this->comm_inter != MPI_COMM_NULL) && (this->comm_inter != MPI_COMM_WORLD)) {

        CHECK(MPI_Comm_free(&this->comm_inter))
        this->comm_intra = MPI_COMM_NULL;
    }
}

//*******************************************************************************
BoundingBox::BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate, const std::shared_ptr<Comm> comm) {

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(dim, 0, comm);

    if (OP_DEBUG)
        opp_printf("Local BoundingBox [provided]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);

}

// For now, only 3D is implemented
//*******************************************************************************
BoundingBox::BoundingBox(const opp_dat node_pos_dat, int dim, const std::shared_ptr<Comm> comm) {
    
    if (dim != 3) {
        opp_abort(std::string("For now, only 3D BoundingBox is implemented"));
    }

    const double* node_pos_data = (const double*)node_pos_dat->data;
    const int node_count = node_pos_dat->set->size + node_pos_dat->set->exec_size + 
                                node_pos_dat->set->nonexec_size;;

    opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
    opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

    // make the bounding box even over the halo regions
    for (int i = 0; i < node_count; i++) {
        minCoordinate.x = std::min(node_pos_data[i * dim + 0], minCoordinate.x);
        minCoordinate.y = std::min(node_pos_data[i * dim + 1], minCoordinate.y);
        minCoordinate.z = std::min(node_pos_data[i * dim + 2], minCoordinate.z);
        maxCoordinate.x = std::max(node_pos_data[i * dim + 0], maxCoordinate.x);
        maxCoordinate.y = std::max(node_pos_data[i * dim + 1], maxCoordinate.y);
        maxCoordinate.z = std::max(node_pos_data[i * dim + 2], maxCoordinate.z);
    }

    if (node_count != 0) {
        minCoordinate.x -= BOUNDING_TOLERENCE;
        minCoordinate.y -= BOUNDING_TOLERENCE;
        minCoordinate.z -= BOUNDING_TOLERENCE;
        maxCoordinate.x += BOUNDING_TOLERENCE;
        maxCoordinate.y += BOUNDING_TOLERENCE;
        maxCoordinate.z += BOUNDING_TOLERENCE;
    }

    this->boundingBox[0] = minCoordinate;
    this->boundingBox[1] = maxCoordinate;

    generateGlobalBoundingBox(dim, node_count, comm);

    if (OP_DEBUG)
        opp_printf("Local BoundingBox [computed]", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->boundingBox[0] .x, this->boundingBox[0] .y, this->boundingBox[0] .z, 
            this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);
}

//*******************************************************************************
void BoundingBox::generateGlobalBoundingBox(int dim, int count, const std::shared_ptr<Comm> comm) {

#ifdef ENABLE_MPI

    const double* localMin = reinterpret_cast<const double*>(&(this->boundingBox[0]));
    double* globalMin = reinterpret_cast<double*>(&(this->globalBoundingBox[0]));    
    MPI_Allreduce(localMin, globalMin, dim, MPI_DOUBLE, MPI_MIN, comm->comm_parent);

    const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
    double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
    MPI_Allreduce(localMax, globalMax, dim, MPI_DOUBLE, MPI_MAX, comm->comm_parent);

    // This is to avoid min and max corrdinates to not have MAX_REAL and MIN_REAL when current rank has no work
    if (count == 0) {
        this->boundingBox[0].x = this->globalBoundingBox[0].x; 
        this->boundingBox[0].y = this->globalBoundingBox[0].y; 
        this->boundingBox[0].z = this->globalBoundingBox[0].z; 

        this->boundingBox[1].x = this->globalBoundingBox[0].x;
        this->boundingBox[1].y = this->globalBoundingBox[0].y;
        this->boundingBox[1].z = this->globalBoundingBox[0].z;
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("Global BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
            this->globalBoundingBox[0].x, this->globalBoundingBox[0].y, this->globalBoundingBox[0].z, 
            this->globalBoundingBox[1].x, this->globalBoundingBox[1].y, this->globalBoundingBox[1].z);
#else
    this->globalBoundingBox = this->boundingBox;
#endif
}

//*******************************************************************************
BoundingBox::~BoundingBox() { };

//*******************************************************************************
const opp_point& BoundingBox::getLocalMin() const {

    return this->boundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getLocalMax() const {

    return this->boundingBox[1];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMin() const {

    return this->globalBoundingBox[0];
}

//*******************************************************************************
const opp_point& BoundingBox::getGlobalMax() const {
    
    return this->globalBoundingBox[1];
}

//*******************************************************************************
bool BoundingBox::isCoordinateInBoundingBox(const opp_point& point) { 

    if (this->boundingBox[0].x > point.x || this->boundingBox[1].x < point.x) 
        return false;
    else if (this->boundingBox[0].y > point.y || this->boundingBox[1].y < point.y) 
        return false;
    else if (this->boundingBox[0].z > point.z || this->boundingBox[1].z < point.z) 
        return false;
    
    return true;
}

//*******************************************************************************
bool BoundingBox::isCoordinateInGlobalBoundingBox(const opp_point& point) { 

    if (this->globalBoundingBox[0].x > point.x || this->globalBoundingBox[1].x < point.x) 
        return false;
    else if (this->globalBoundingBox[0].y > point.y || this->globalBoundingBox[1].y < point.y) 
        return false;
    else if (this->globalBoundingBox[0].z > point.z || this->globalBoundingBox[1].z < point.z) 
        return false;
    
    return true;
}


//*******************************************************************************
// This will contain mappings for halo indices too
GlobalToLocalCellIndexMapper::GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat) {

    // if (OP_DEBUG) 
    if (OPP_rank == 0)            
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping start");
        
    globalToLocalCellIndexMap.clear();

    const opp_set cells_set = global_cell_id_dat->set;
    int size_inc_halo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;

    for (int i = 0; i < size_inc_halo; i++) {

        int glbIndex = ((int*)global_cell_id_dat->data)[i];
        
        globalToLocalCellIndexMap.insert(std::make_pair(glbIndex, i));
    }
    
    // if (OP_DEBUG) 
    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "generateGlobalToLocalCellIndexMapping end");
}

//*******************************************************************************
GlobalToLocalCellIndexMapper::~GlobalToLocalCellIndexMapper() {
    
    globalToLocalCellIndexMap.clear();

    if (OPP_rank == 0)
        opp_printf("GlobalToLocalCellIndexMapper", "destroyGlobalToLocalCellIndexMapping");
}

//*******************************************************************************
int GlobalToLocalCellIndexMapper::map(const int globalIndex) {

    if (globalIndex == MAX_CELL_INDEX)
        return MAX_CELL_INDEX;

    auto it = globalToLocalCellIndexMap.find(globalIndex);
    if (it != globalToLocalCellIndexMap.end())
        return it->second;
    
    opp_printf("GlobalToLocalCellIndexMapper", "Error... local cell index not found for global index %d", globalIndex);
    return MAX_INT;
}


//*******************************************************************************
CellMapper::CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, const std::shared_ptr<Comm> comm) 
    : boundingBox(boundingBox), gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing), 
        minGlbCoordinate(boundingBox->getGlobalMin()), comm(comm) {
    
    const opp_point& minGblCoordinate = boundingBox->getGlobalMin();
    const opp_point& maxGblCoordinate = boundingBox->getGlobalMax();

    int ax = 0, ay = 0, az = 0;

    { // removed this and added below due to decimal point issues
        // this->globalGridDims.x = std::ceil((maxGblCoordinate.x - minGblCoordinate.x) * oneOverGridSpacing);
        // this->globalGridDims.y = std::ceil((maxGblCoordinate.y - minGblCoordinate.y) * oneOverGridSpacing);
        // this->globalGridDims.z = std::ceil((maxGblCoordinate.z - minGblCoordinate.z) * oneOverGridSpacing);   
    }
    {
        for (double z = minGblCoordinate.z; z < maxGblCoordinate.z; z += this->gridSpacing) az++;
        for (double y = minGblCoordinate.y; y < maxGblCoordinate.y; y += this->gridSpacing) ay++;
        for (double x = minGblCoordinate.x; x < maxGblCoordinate.x; x += this->gridSpacing) ax++; 
        this->globalGridDims.x = ax + 1;
        this->globalGridDims.y = ay + 1;
        this->globalGridDims.z = az + 1; 
    }

    if (OPP_rank == OPP_ROOT)
        opp_printf("CellMapper", "Global Grid Size - [%d %d %d] gridSpacing [%2.10lE]", 
            this->globalGridDims.x, this->globalGridDims.y, this->globalGridDims.z, this->gridSpacing); 
    
    const opp_point& minLocalCoordinate = boundingBox->getLocalMin();
    const opp_point& maxLocalCoordinate = boundingBox->getLocalMax();

    // Find the local ranks grid start indices
    ax = 0; ay = 0; az = 0;
    for (double z = minGblCoordinate.z; (z < minLocalCoordinate.z); z += this->gridSpacing) az++;
    for (double y = minGblCoordinate.y; (y < minLocalCoordinate.y); y += this->gridSpacing) ay++; 
    for (double x = minGblCoordinate.x; (x < minLocalCoordinate.x); x += this->gridSpacing) ax++; 
    this->localGridStart.x = ax == 0 ? ax : (ax - 1);
    this->localGridStart.y = ay == 0 ? ay : (ay - 1);
    this->localGridStart.z = az == 0 ? az : (az - 1);         

    // Find the local ranks grid end indices
    ax = 0; ay = 0; az = 0;
    for (double z = minGblCoordinate.z; (z <= maxLocalCoordinate.z); z += this->gridSpacing) az++; 
    for (double y = minGblCoordinate.y; (y <= maxLocalCoordinate.y); y += this->gridSpacing) ay++; 
    for (double x = minGblCoordinate.x; (x <= maxLocalCoordinate.x); x += this->gridSpacing) ax++; 
    this->localGridEnd.x = this->globalGridDims.x == ax ? ax : (ax + 1);
    this->localGridEnd.y = this->globalGridDims.y == ay ? ay : (ay + 1);
    this->localGridEnd.z = this->globalGridDims.z == az ? az : (az + 1);     

    if (OP_DEBUG)
        opp_printf("CellMapper", "Local Grid - Min[%d %d %d] Max[%d %d %d]", 
            this->localGridStart.x, this->localGridStart.y, this->localGridStart.z, 
            this->localGridEnd.x, this->localGridEnd.y, this->localGridEnd.z); 
}

//*******************************************************************************
CellMapper::~CellMapper() { 

#ifdef ENABLE_MPI
    CHECK(MPI_Win_free(&this->win_structMeshToCellMapping))
    this->structMeshToCellMapping = nullptr;

    CHECK(MPI_Win_free(&this->win_structMeshToRankMapping))
    this->structMeshToRankMapping = nullptr;           
#else
    delete[] this->structMeshToCellMapping;
#endif
};

//*******************************************************************************
opp_point CellMapper::getCentroidOfBox(const opp_point& coordinate) { 

    opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);
    const opp_point& maxCoordinate = boundingBox->getGlobalMax();

    constexpr int DIM = 3; // Hardcoded for now
    switch (DIM) {
        case 1:
            ASSIGN_CENTROID_TO_DIM(x); 
            break;
        case 2:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            break;
        case 3:
            ASSIGN_CENTROID_TO_DIM(x);
            ASSIGN_CENTROID_TO_DIM(y);
            ASSIGN_CENTROID_TO_DIM(z);
            break;
        default:
            std::cerr << "Error getCentroidOfBox: Dimension invalid " << DIM << std::endl;
    }

    return centroid;
}

//*******************************************************************************
// Returns the global cell index
size_t CellMapper::findStructuredCellIndex(const opp_point& position) { 

    // Perform the calculations in higher precision (double)
    double xDiff = position.x - this->minGlbCoordinate.x;
    double yDiff = position.y - this->minGlbCoordinate.y;
    double zDiff = position.z - this->minGlbCoordinate.z;

    xDiff = xDiff * this->oneOverGridSpacing;
    yDiff = yDiff * this->oneOverGridSpacing;
    zDiff = zDiff * this->oneOverGridSpacing;

    // Round to the nearest integer to minimize rounding errors
    const int xIndex = static_cast<int>(xDiff);
    const int yIndex = static_cast<int>(yDiff);
    const int zIndex = static_cast<int>(zDiff);

    // Calculate the cell index mapping index
    size_t index = ((size_t)(xIndex) + (size_t)(yIndex * globalGridDims.x) + 
                (size_t)(zIndex * globalGridDims.x * globalGridDims.y));

    if (index >= globalGridSize) {
        // opp_printf("CellMapper", "Error index %zu generated is larger than globalGridSize %zu", 
        //     index, globalGridSize);
        return MAX_CELL_INDEX;
    }

    return index;
}

//*******************************************************************************
// Returns the global cell index
int CellMapper::findClosestCellIndex(const size_t& structCellIdx) { 
    
    if (OP_DEBUG) 
    {
        if (structCellIdx >= globalGridSize) {
            opp_printf("findClosestCellIndex", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                structCellIdx, globalGridSize);
            return MAX_CELL_INDEX;
        }
    }
    
    return this->structMeshToCellMapping[structCellIdx];
}

//*******************************************************************************
// Returns the rank of the cell
int CellMapper::findClosestCellRank(const size_t& structCellIdx) { 

#ifdef ENABLE_MPI 
    if (OP_DEBUG) 
    {
        if (structCellIdx >= globalGridSize) {
            opp_printf("findClosestCellRank", "Warning returning MAX - structCellIdx=%zu globalGridSize=%zu",
                structCellIdx, globalGridSize);
            return MAX_CELL_INDEX;
        }
    }

    return this->structMeshToRankMapping[structCellIdx];
#else
    return OPP_rank;
#endif
}

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) {

#ifdef ENABLE_MPI
    waitBarrier();

    if (comm->rank_intra == 0) { // comm->size_intra > 1 && 

        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));

        if (OP_DEBUG) 
            opp_printf("CellMapper", "reduceInterNodeMappings Start size_inter %d", comm->size_inter);

        std::string sendLog = "Send Counts: ", recvLog = "Recv Counts: ";

        int totalSendCount = 0, totalRecvCount = 0;
        std::vector<int> sendCounts(comm->size_inter);
        std::vector<int> sendDisplacements(comm->size_inter);
        std::vector<int> recvCounts(comm->size_inter);
        std::vector<int> recvDisplacements(comm->size_inter);

        for (int i = 0; i < comm->size_inter; ++i) {
            const int remainder = (int)(this->globalGridSize % comm->size_inter);
            sendCounts[i] = (this->globalGridSize / comm->size_inter) + (i < remainder ? 1 : 0);
            sendDisplacements[i] = totalSendCount;
            totalSendCount += sendCounts[i];

            if (OP_DEBUG) 
                sendLog += std::to_string(sendCounts[i]) + "|Disp:" + std::to_string(sendDisplacements[i]) + " ";
        }

        for (int i = 0; i < comm->size_inter; ++i) {
            recvCounts[i] = sendCounts[comm->rank_inter];
            recvDisplacements[i] = totalRecvCount;
            totalRecvCount += recvCounts[i];

            if (OP_DEBUG) 
                recvLog += std::to_string(recvCounts[i]) + "|Disp:" + std::to_string(recvDisplacements[i]) + " ";
        }

        if (OP_DEBUG) 
        {
            opp_printf("CellMapper", "reduceInterNodeMappings totalSendCount %d : %s", totalSendCount, sendLog.c_str());
            opp_printf("CellMapper", "reduceInterNodeMappings totalRecvCount %d : %s", totalRecvCount, recvLog.c_str());
        }

        std::vector<int> cellMappingsRecv(totalRecvCount, MAX_CELL_INDEX-1);
        std::vector<int> ranksRecv(totalRecvCount, MAX_CELL_INDEX-1);

        const int INT_COUNT_PER_MSG = 10000;

        std::vector<MPI_Request> sendRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalSendCount = sendCounts[rank];
            int sendCount = 0, alreadySentCount = 0;

            int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {
                
                if (i != blocks-1) 
                    sendCount = INT_COUNT_PER_MSG;
                else 
                    sendCount = totalSendCount - alreadySentCount;

                // opp_printf("SEND", "blocks %d|%d disp %d sendCounts %d to rank %d | %d %d %d %d %d %d %d %d %d %d", 
                //     i, blocks,
                //     sendDisplacements[rank] + i * INT_COUNT_PER_MSG, sendCount, inter_ranks[rank],
                //     sb[0], sb[1], sb[2], sb[3], sb[4], sb[5], sb[6], sb[7], sb[8], sb[9]);

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                        rank, (10000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, 
                        rank, (20000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                alreadySentCount += sendCount;
            }
        }

        std::vector<MPI_Request> recvRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalRecvCount = recvCounts[rank];
            int recvCount = 0, alreadyRecvCount = 0;

            int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {

                if (i != blocks-1) 
                    recvCount = INT_COUNT_PER_MSG;
                else 
                    recvCount = totalRecvCount - alreadyRecvCount;

                // opp_printf("RECV", "blocks %d|%d disp %d recvCounts %d from rank %d", i, blocks,
                //     recvDisplacements[rank] + i * INT_COUNT_PER_MSG, recvCount, inter_ranks[rank]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(ranksRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                        rank, (10000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(cellMappingsRecv[recvDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount, MPI_INT, 
                        rank, (20000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                alreadyRecvCount += recvCount;
            }
        }

        // opp_printf("DISPLACEMENTS", "send %d recv %d", sendDisplacements[comm->rank_inter], recvDisplacements[comm->rank_inter]);

        // Copy own data
        for (int i = 0; i < recvCounts[comm->rank_inter]; i++) {

            ranksRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i];
            cellMappingsRecv[recvDisplacements[comm->rank_inter] + i] = this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i];
        }

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);

        // printStructuredMesh(std::string("RECV_MAPPING") + std::to_string(callID), cellMappingsRecv.data(), totalRecvCount);
        // printStructuredMesh(std::string("RECV_RANKS") + std::to_string(callID), ranksRecv.data(), totalRecvCount);
        // MPI_Barrier(comm->comm_inter);

        const size_t recvCount = recvCounts[0];
        for (size_t i = 0; i < recvCount; i++) { // reduce to get common mapping on all inter node ranks

            int cellIndex = MAX_CELL_INDEX;
            for (int r = 0; r < comm->size_inter; r++) {

                if (cellIndex > cellMappingsRecv[i + r * recvCount]) {

                    cellIndex = cellMappingsRecv[i + r * recvCount];
                    cellMappingsRecv[i] = cellIndex;
                    ranksRecv[i] = ranksRecv[i + r * recvCount];
                }
            }
        }

        MPI_Barrier(comm->comm_inter);

        // printStructuredMesh(std::string("MAPPING_COMP") + std::to_string(callID), ranksRecv.data(), recvCount);
        // MPI_Barrier(comm->comm_inter);

        sendRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalSendCount = recvCount;
            int sendCount = 0, alreadySentCount = 0;

            int blocks = (totalSendCount / INT_COUNT_PER_MSG) + ((totalSendCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {
                
                if (i != blocks-1) 
                    sendCount = INT_COUNT_PER_MSG;
                else 
                    sendCount = totalSendCount - alreadySentCount;

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(cellMappingsRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (30000 + i), 
                        comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                sendRequests.emplace_back(MPI_Request());
                MPI_Isend(&(ranksRecv[i * INT_COUNT_PER_MSG]), sendCount, MPI_INT, rank, (40000 + i), 
                        comm->comm_inter, &sendRequests[sendRequests.size() - 1]); 

                alreadySentCount += sendCount;
            }
        }

        // Since we are about to write to data mapped by windows, release the lock here
        CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
        CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));

        recvRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) 
        {
            if (rank == comm->rank_inter) continue;

            int totalRecvCount = sendCounts[rank];
            int recvCount2 = 0, alreadyRecvCount = 0;

            int blocks = (totalRecvCount / INT_COUNT_PER_MSG) + ((totalRecvCount % INT_COUNT_PER_MSG == 0) ? 0 : 1);

            for (int i = 0; i < blocks; i++) {

                if (i != blocks-1) 
                    recvCount2 = INT_COUNT_PER_MSG;
                else 
                    recvCount2 = totalRecvCount - alreadyRecvCount;

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(this->structMeshToCellMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                        rank, (30000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                recvRequests.emplace_back(MPI_Request());
                MPI_Irecv(&(this->structMeshToRankMapping[sendDisplacements[rank] + i * INT_COUNT_PER_MSG]), recvCount2, MPI_INT, 
                        rank, (40000 + i), comm->comm_inter, &recvRequests[recvRequests.size() - 1]);

                alreadyRecvCount += recvCount2;
            }
        }

        // Copy own data
        for (int i = 0; i < sendCounts[comm->rank_inter]; i++) {

            this->structMeshToRankMapping[sendDisplacements[comm->rank_inter] + i] = ranksRecv[i];
            this->structMeshToCellMapping[sendDisplacements[comm->rank_inter] + i] = cellMappingsRecv[i];
        }

        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE);
        MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE);
        MPI_Barrier(comm->comm_inter);
        
        if (OP_DEBUG) 
            opp_printf("CellMapper", "reduceInterNodeMappings END");
    }

    waitBarrier();
#endif
}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OP_DEBUG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat);

#ifdef ENABLE_MPI
    if (comm->rank_intra == 0) {
        
        CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));

        for (size_t i = 0; i < globalGridSize; i++) {
            
            if (this->structMeshToCellMapping[i] != MAX_CELL_INDEX)
                this->structMeshToCellMapping[i] = (-1 * this->structMeshToCellMapping[i]);
        }

        CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
    }

    MPI_Win_fence(0, this->win_structMeshToCellMapping); 
#endif

    for (size_t i = 0; i < globalGridSize; i++) {
        
        if (this->structMeshToRankMapping[i] == OPP_rank) {
            
            const int globalCID = (-1 * this->structMeshToCellMapping[i]);
            
            if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {
                
                const int localCID = globalToLocalCellIndexMapper.map(globalCID);
                
                if (localCID != MAX_CELL_INDEX) {
                    this->structMeshToCellMapping[i] = localCID;
                }
                else {
                    opp_printf("CellMapper", "Error at convertToLocalMappings : structMeshToCellMapping at index %d is invalid [gcid:%d]", 
                        i, this->structMeshToCellMapping[i]);
                }
            }
        }
    }

#ifdef ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);

    if (comm->rank_intra == 0) {

        MPI_Allreduce(MPI_IN_PLACE, this->structMeshToCellMapping, globalGridSize, MPI_INT, MPI_MAX, comm->comm_inter);
    }

    waitBarrier();
#endif
    if (OP_DEBUG) 
    // if (OPP_rank == 0)
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::enrichStructuredMesh(const int index, const int cell_index, const int rank) {

#ifdef ENABLE_MPI
    MPI_Put(&cell_index, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToCellMapping);
    MPI_Put(&rank, 1, MPI_INT, 0, index, 1, MPI_INT, this->win_structMeshToRankMapping);
#else
    this->structMeshToCellMapping[index] = cell_index;
#endif
}

//*******************************************************************************
void CellMapper::printStructuredMesh(const std::string msg, int *array, size_t size, bool printToFile) {
    
    // if (!OP_DEBUG)
    //     return;

    if (!printToFile) 
    {
        opp_printf("structMeshToCellMapping", "%s - size=%zu", msg.c_str(), size);

        for (size_t i = 0; i < size; i++) {
            printf("%zu|%d ", i, array[i]);
            if (i % 50 == 0) printf("\n");
        }   
        printf("\n");
    }
    else
    {
        const std::string file_name = std::string("files/struct_com") + std::to_string(OPP_comm_size) + "_r" + 
                                        std::to_string(OPP_rank) + "_" + msg; 

        std::ofstream outFile(file_name);

        if (!outFile.is_open()) {
            opp_printf("printStructuredMesh", "can't open file %s\n", file_name.c_str());
            opp_abort("printStructuredMesh - can't open file");
        }

        for (size_t i = 0; i < size; i++) {

            outFile << i << "|";

            int value = array[i];
            if (value != MAX_CELL_INDEX) {
                outFile << value << " ";
            }
            else {
                outFile << "X ";
            }

            if ((i+1) % 50 == 0) 
                outFile << "\n";
        }   
        outFile << "\n";
        outFile.close();
    }
}

//*******************************************************************************
void CellMapper::createStructMeshMappingArrays() {

    globalGridSize = (size_t)(globalGridDims.x * globalGridDims.y * globalGridDims.z);

#ifdef ENABLE_MPI

    MPI_Barrier(MPI_COMM_WORLD);

    // One per shared memory (node)
    
    // create CELL INDEX array
    {
        const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

        CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                    (void *)&this->structMeshToCellMapping, &this->win_structMeshToCellMapping))

        MPI_Aint allocatedSize = 0;
        int disp = 0;
        CHECK(MPI_Win_shared_query(this->win_structMeshToCellMapping, 0, &allocatedSize, &disp,
                                    (void *)&this->structMeshToCellMapping))

        if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
            opp_abort(std::string("Pointer to incorrect size in MPI structMeshToCellMapping"));
        }
        if (disp != sizeof(int)) {
            opp_abort(std::string("Invalid displacement unit in MPI structMeshToCellMapping"));
        }

        if (comm->rank_intra == 0) {
            for (size_t i = 0; i < globalGridSize; i++)
                this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
        }

        MPI_Win_fence(0, this->win_structMeshToCellMapping);
    }

    // create RANK array
    {
        const MPI_Aint size = (comm->rank_intra == 0) ? globalGridSize * sizeof(int) : 0;

        CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                    (void *)&this->structMeshToRankMapping, &this->win_structMeshToRankMapping))

        MPI_Aint allocatedSize = 0;
        int disp = 0;
        CHECK(MPI_Win_shared_query(this->win_structMeshToRankMapping, 0, &allocatedSize, &disp,
                                    (void *)&this->structMeshToRankMapping))

        if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
            opp_abort(std::string("Pointer to incorrect size in MPI structMeshToRankMapping"));
        }
        if (disp != sizeof(int)) {
            opp_abort(std::string("Invalid displacement unit in MPI structMeshToRankMapping"));
        }

        if (comm->rank_intra == 0) {
            for (size_t i = 0; i < globalGridSize; i++)
                this->structMeshToRankMapping[i] = MAX_CELL_INDEX;
        }

        MPI_Win_fence(0, this->win_structMeshToRankMapping);
    }
#else
    this->structMeshToCellMapping = new int[globalGridSize];
    for (size_t i = 0; i < globalGridSize; i++)
        this->structMeshToCellMapping[i] = MAX_CELL_INDEX;
#endif
}

//*******************************************************************************
void CellMapper::waitBarrier() {

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_fence(0, this->win_structMeshToCellMapping); 
    MPI_Win_fence(0, this->win_structMeshToRankMapping); 
}

//*******************************************************************************
void CellMapper::lockWindows() {

    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, this->win_structMeshToRankMapping));
}


void CellMapper::unlockWindows() {

    CHECK(MPI_Win_unlock(0, this->win_structMeshToCellMapping));
    CHECK(MPI_Win_unlock(0, this->win_structMeshToRankMapping));
}
