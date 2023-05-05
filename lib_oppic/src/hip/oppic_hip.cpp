
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

#include <oppic_hip.h>

//****************************************
void oppic_init(int argc, char **argv, opp::Params* params)
{
    oppic_init_core(argc, argv, params);
    cutilDeviceInit(argc, argv);

    // cutilSafeCall(hipDeviceSetCacheConfig(hipFuncCachePreferShared));
    // cutilSafeCall(hipDeviceSetSharedMemConfig(hipSharedMemBankSizeEightByte));

    OP_auto_soa = 1; // TODO : Make this configurable with args
    OP_auto_sort = 1;
}

//****************************************
void oppic_exit()
{
    oppic_cuda_exit();
    oppic_exit_core();
}

//****************************************
void oppic_cuda_exit() 
{
    if (!OP_hybrid_gpu)
        return;

    for (auto& a : oppic_maps) 
    {
        cutilSafeCall(hipFree(a->map_d));
    }

    for (auto& a : oppic_dats) 
    {
        // cutilSafeCall(hipFree(a->data_d));
        if (a->thrust_int) delete a->thrust_int;
        if (a->thrust_real) delete a->thrust_real;
        if (a->thrust_int_sort) delete a->thrust_int_sort;
        if (a->thrust_real_sort) delete a->thrust_real_sort;
    }

}

//****************************************
oppic_set oppic_decl_set(int size, char const *name)
{
    return oppic_decl_set_core(size, name);
}

//****************************************
oppic_map oppic_decl_map(oppic_set from, oppic_set to, int dim, int *imap, char const *name)
{
    oppic_map map = oppic_decl_map_core(from, to, dim, imap, name);

    int set_size = map->from->size + map->from->exec_size;
    int *temp_map = (int *)malloc(map->dim * set_size * sizeof(int));

    for (int i = 0; i < map->dim; i++) 
    {
        for (int j = 0; j < set_size; j++) 
        {
            temp_map[i * set_size + j] = map->map[map->dim * j + i];
        }
    }

    int copy_size = map->dim * set_size * sizeof(int);
    oppic_cpHostToDevice((void **)&(map->map_d), (void **)&(temp_map), copy_size, copy_size, true);
    free(temp_map);

    return map;
}

//****************************************
oppic_dat oppic_decl_dat(oppic_set set, int dim, opp_data_type dtype, char *data, char const *name)
{
    std::string type = "";
    int size = -1;

    getDatTypeSize(dtype, type, size);
    
    oppic_dat dat = oppic_decl_dat_core(set, dim, type.c_str(), size, data, name);

    oppic_create_device_arrays(dat);

    oppic_upload_dat(dat);

    return dat;
}

//****************************************
oppic_map oppic_decl_map_txt(oppic_set from, oppic_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)oppic_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    oppic_map map = oppic_decl_map(from, to, dim, map_data, name);

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
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, mapping);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, idx, map, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_dat dat, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(dat, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_map data_map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, acc, mapping);
}
oppic_arg oppic_arg_dat(oppic_map data_map, int idx, oppic_map map, oppic_access acc, opp_mapping mapping)
{
    return oppic_arg_dat_core(data_map, idx, map, acc, mapping);
}

//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg oppic_arg_gbl(int *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}
oppic_arg oppic_arg_gbl(const bool *data, int dim, char const *typ, oppic_access acc)
{
    return oppic_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
oppic_set oppic_decl_particle_set(int size, char const *name, oppic_set cells_set)
{
    oppic_set set = oppic_decl_particle_set_core(size, name, cells_set);

    cutilSafeCall(hipMalloc((void**)&(set->particle_remove_count_d), sizeof(int)));
    cutilSafeCall(hipDeviceSynchronize());
  
    return set;
}
oppic_set oppic_decl_particle_set(char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set(0, name, cells_set);
}

//****************************************
oppic_dat oppic_decl_particle_dat(oppic_set set, int dim, opp_data_type dtype, char *data, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    oppic_dat dat = oppic_decl_particle_dat_core(set, dim, type.c_str(), size, data, name, cell_index);

    oppic_create_device_arrays(dat);

    oppic_upload_dat(dat);

    return dat;
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
{ TRACE_ME;

    bool need_resizing = (particles_set->set_capacity < (particles_set->size + num_particles_to_insert)) ? true : false;

    if (OP_DEBUG) printf("\toppic_increase_particle_count need_resizing %s\n", need_resizing ? "YES" : "NO");

    if (need_resizing)
        oppic_download_particle_set(particles_set); // TODO : We should be able to do a device to device copy instead of getting to host

    oppic_increase_particle_count_core(particles_set, num_particles_to_insert);

    if (need_resizing)
    {
        for (oppic_dat& current_dat : *(particles_set->particle_dats))
        {
            if (OP_DEBUG) printf("\toppic_increase_particle_count | dat [%s]\n", current_dat->name);

            // TODO : We might be able to copy only the old data from device to device!

            oppic_create_device_arrays(current_dat, true);

            oppic_upload_dat(current_dat);

            current_dat->dirty_hd = Dirty::NotDirty;
        }        
    }

}

//****************************************
void oppic_init_particle_move(oppic_set set)
{ TRACE_ME;

    oppic_init_particle_move_core(set);

    cutilSafeCall(hipMemcpy(set->particle_remove_count_d, &(set->particle_remove_count), sizeof(int), hipMemcpyHostToDevice));
}

//****************************************
bool oppic_finalize_particle_move(oppic_set set)
{ TRACE_ME;

    hipMemcpy(&(set->particle_remove_count), set->particle_remove_count_d, sizeof(int), hipMemcpyDeviceToHost);

    if (OP_DEBUG) printf("\toppic_finalize_particle_move set [%s] with particle_remove_count [%d]\n", set->name, set->particle_remove_count);
    
    // oppic_finalize_particle_move_cuda(set); // This makes device-host-device copies, auto sorting takes less time!
    
    oppic_finalize_particle_move_core(set); // OPP_auto_sorting should be true

    if (OP_auto_sort == 1)
    {
        if (OP_DEBUG) printf("\toppic_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        oppic_particle_sort(set);
    }

    return true;
}

//****************************************
void particle_sort_cuda(oppic_set set);

//****************************************
void oppic_particle_sort(oppic_set set)
{ TRACE_ME;

    particle_sort_cuda(set);
}

//****************************************
void oppic_print_dat_to_txtfile(oppic_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    if (OP_DEBUG) printf("\toppic_print_dat_to_txtfile writing file [%s]\n", file_name_suffix);

    if (dat->dirty_hd == Dirty::Host) oppic_download_dat(dat);

    std::string prefix = std::string(file_name_prefix) + "_c";
    oppic_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void oppic_print_map_to_txtfile(oppic_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    oppic_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

//****************************************
void oppic_dump_dat(oppic_dat dat)
{
    if (dat->dirty_hd == Dirty::Host) oppic_download_dat(dat);

    oppic_dump_dat_core(dat);
}

//****************************************
// DEVICE->HOST | this invalidates what is in the HOST
void oppic_download_dat(oppic_dat dat) 
{ TRACE_ME;

    size_t set_size = dat->set->set_capacity;
    if (strstr(dat->type, ":soa") != NULL || (OP_auto_soa && dat->dim > 1)) 
    {
        if (OP_DEBUG) printf("\toppic_download_dat GPU->CPU SOA | %s\n", dat->name);

        char *temp_data = (char *)malloc(dat->size * set_size * sizeof(char));
        cutilSafeCall(hipMemcpy(temp_data, dat->data_d, set_size * dat->size, hipMemcpyDeviceToHost));
        
        int element_size = dat->size / dat->dim;
        for (int i = 0; i < dat->dim; i++) 
        {
            for (int j = 0; j < set_size; j++) 
            {
                for (int c = 0; c < element_size; c++) 
                {
                    dat->data[dat->size * j + element_size * i + c] = temp_data[element_size * i * set_size + element_size * j + c];
                }
            }
        }
        free(temp_data);
    } 
    else 
    {
        if (OP_DEBUG) printf("\toppic_download_dat GPU->CPU NON-SOA| %s\n", dat->name);

        cutilSafeCall(hipMemcpy(dat->data, dat->data_d, set_size * dat->size, hipMemcpyDeviceToHost));
    }

    dat->dirty_hd = Dirty::NotDirty;
}

//****************************************
// HOST->DEVICE | this invalidates what is in the DEVICE
void oppic_upload_dat(oppic_dat dat)
{ TRACE_ME;

    int set_capacity = dat->set->set_capacity;

    if (strstr(dat->type, ":soa") != NULL || (OP_auto_soa && dat->dim > 1)) 
    {
        if (OP_DEBUG) printf("\toppic_upload_dat CPU->GPU SOA | %s\n", dat->name);

        char *temp_data = (char *)malloc(dat->size * set_capacity * sizeof(char));
        int element_size = dat->size / dat->dim;

        for (int i = 0; i < dat->dim; i++) 
        {
            for (int j = 0; j < set_capacity; j++) 
            {
                for (int c = 0; c < element_size; c++) 
                {
                    temp_data[element_size * i * set_capacity + element_size * j + c] = dat->data[dat->size * j + element_size * i + c];
                }
            }
        }

        oppic_cpHostToDevice((void **)&(dat->data_d), (void **)&(temp_data), (dat->size * set_capacity), (dat->size * set_capacity), false);
        free(temp_data);
    } 
    else 
    {
        if (OP_DEBUG) printf("\toppic_upload_dat CPU->GPU NON-SOA| %s\n", dat->name);
        oppic_cpHostToDevice((void **)&(dat->data_d), (void **)&(dat->data), (dat->size * set_capacity), (dat->size * set_capacity), false);
    }

    dat->dirty_hd = Dirty::NotDirty;
}

//****************************************
// DEVICE -> HOST
void oppic_download_particle_set(oppic_set particles_set)
{

    if (OP_DEBUG) printf("\toppic_download_particle_set set [%s]\n", particles_set->name);

    for (oppic_dat& current_dat : *(particles_set->particle_dats))
    {
        if (current_dat->data_d == NULL)
        {
            if (OP_DEBUG) printf("\toppic_download_particle_set device pointer is NULL in dat [%s]\n", current_dat->name);
            continue;
        }
        if (current_dat->dirty_hd != Dirty::Host)
        {
            if (OP_DEBUG) printf("\toppic_download_particle_set host is not dirty in dat [%s]\n", current_dat->name);
            continue;
        }
        oppic_download_dat(current_dat);
    }  
}

//****************************************
// HOST->DEVICE
void oppic_upload_particle_set(oppic_set particles_set, bool realloc)
{ 

    if (OP_DEBUG) printf("\toppic_upload_particle_set set [%s]\n", particles_set->name);

    for (oppic_dat& current_dat : *(particles_set->particle_dats))
    {
        if (realloc)
        {
            // TODO : CONVERT TO THRUST VECTORS, WILL BREAK IF NOT 
            if (current_dat->data_d != NULL) 
                cutilSafeCall(hipFree(current_dat->data_d));

            cutilSafeCall(hipMalloc(&(current_dat->data_d), particles_set->set_capacity * current_dat->size));
        }

        oppic_upload_dat(current_dat);
    }  
}

//****************************************
int oppic_mpi_halo_exchanges_cuda(oppic_set set, int nargs, oppic_arg *args);
void oppic_mpi_set_dirtybit_cuda(int nargs, oppic_arg *args);

int oppic_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device)
{
    return (device == Device_CPU) ? oppic_mpi_halo_exchanges(set, nargs, args) : oppic_mpi_halo_exchanges_cuda(set, nargs, args);
}

//****************************************
// DEVICE->HOST copy of Dirty::Host dats 
int oppic_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args) 
{ TRACE_ME;
    for (int n = 0; n < nargs; n++)
    {
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == Dirty::Host) 
        {
            oppic_download_dat(args[n].dat);
            args[n].dat->dirty_hd = Dirty::NotDirty;
        }
    }
    return set->size;
}

//****************************************
// HOST->DEVICE copy of Dirty::Device dats
int oppic_mpi_halo_exchanges_cuda(oppic_set set, int nargs, oppic_arg *args) 
{ TRACE_ME;
    for (int n = 0; n < nargs; n++)
    { 
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == Dirty::Device) 
        { 
            oppic_upload_dat(args[n].dat);
            args[n].dat->dirty_hd = Dirty::NotDirty;
        }
    }
    return set->size;
}

//****************************************
void oppic_mpi_set_dirtybit_grouped(int nargs, oppic_arg *args, DeviceType device)
{
    return (device == Device_CPU) ? oppic_mpi_set_dirtybit(nargs, args) : oppic_mpi_set_dirtybit_cuda(nargs, args);
}

//****************************************
// Set Dirty::Device, making HOST data to be clean
void oppic_mpi_set_dirtybit(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
            (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
            args[n].acc == OP_RW)) 
        {
            if (OP_DEBUG) printf("\toppic_mpi_set_dirtybit Setting Dirty::Device| %s\n", args[n].dat->name);
            args[n].dat->dirty_hd = Dirty::Device;
        }
    }
}

//****************************************
// Set Dirty::Host, making DEVICE data to be clean
void oppic_mpi_set_dirtybit_cuda(int nargs, oppic_arg *args) 
{
    for (int n = 0; n < nargs; n++) 
    {
        if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
            (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
            args[n].acc == OP_RW)) 
        {
            if (OP_DEBUG) printf("\top_mpi_set_dirtybit_cuda Setting Dirty::Host| %s\n", args[n].dat->name);
            args[n].dat->dirty_hd = Dirty::Host;
        }
    }
}

//****************************************
// Set the complete dat to zero (upto array capacity)
void oppic_reset_dat(oppic_dat dat, char* val, opp_reset reset)
{
    UNUSED(val);

    if (OP_DEBUG) printf("\toppic_reset_dat dat [%s]\n", dat->name);

    int set_size = dat->set->set_capacity;

    cutilSafeCall(hipMemset((double*)(dat->data_d), 0, dat->size * set_size));

    dat->dirty_hd = Dirty::Host;
}







// **************************************** UTILITY FUNCTIONS ****************************************

void __cudaSafeCall(hipError_t err, const char *file, const int line) 
{
    if (hipSuccess != err) 
    {
        fprintf(stderr, "%s(%i) : cutilSafeCall() Runtime API error : %s.\n", file, line, hipGetErrorString(err));
        exit(-1);
    }
}

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line) 
{
    hipError_t err = hipGetLastError();
    if (hipSuccess != err) 
    {
        fprintf(stderr, "%s(%i) : cutilCheckMsg() error : %s : %s.\n", file, line, errorMessage, hipGetErrorString(err));
        exit(-1);
    }
}


void oppic_cpHostToDevice(void **data_d, void **data_h, int copy_size, int alloc_size, bool create_new) 
{
    if (create_new)
    {
        if (*data_d != NULL) cutilSafeCall(hipFree(*data_d));
        cutilSafeCall(hipMalloc(data_d, alloc_size));
    }

    cutilSafeCall(hipMemcpy(*data_d, *data_h, copy_size, hipMemcpyHostToDevice));
    cutilSafeCall(hipDeviceSynchronize());
}

void oppic_create_device_arrays(oppic_dat dat, bool create_new)
{
    if (OP_DEBUG) printf("oppic_create_device_arrays %s %s\n", dat->name, dat->type);

    if (strcmp(dat->type, "double") == 0)
    {
        if (dat->set->size > 0)
        {
            if (create_new && dat->thrust_real)
            {
                delete dat->thrust_real;
                delete dat->thrust_real_sort;
            }

            dat->thrust_real = new thrust::device_vector<double>(dat->set->set_capacity * dat->dim);
            dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_real->data());

            dat->thrust_real_sort = new thrust::device_vector<double>(dat->set->set_capacity * dat->dim);
        } 
    } 
    else if (strcmp(dat->type, "int") == 0 )
    {
        if (dat->set->size > 0)
        {
            if (create_new && dat->thrust_int)
            {
                delete dat->thrust_int;
                delete dat->thrust_int_sort;
            }

            dat->thrust_int = new thrust::device_vector<int>(dat->set->set_capacity * dat->dim);
            dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_int->data());

            dat->thrust_int_sort = new thrust::device_vector<int>(dat->set->set_capacity * dat->dim);
        } 
    }
    else
    {
        std::cerr << "oppic_decl_dat CUDA not implemented for type: " << dat->type << " dat name: " << dat->name << std::endl;
        exit(-1);
    }
}

void cutilDeviceInit(int argc, char **argv) 
{
    (void)argc;
    (void)argv;
    int deviceCount;

    cutilSafeCall(hipGetDeviceCount(&deviceCount));
    if (deviceCount == 0) 
    {
        printf("cutil error: no devices supporting CUDA\n");
        exit(-1);
    }

    // Test we have access to a device
    float *test;
    hipError_t err = hipMalloc((void **)&test, sizeof(float));
    if (err != hipSuccess) OP_hybrid_gpu = 0;
    else  OP_hybrid_gpu = 1;

    if (OP_hybrid_gpu) 
    {
        hipFree(test);

        // cutilSafeCall(hipDeviceSetCacheConfig(hipFuncCachePreferL1));

        int deviceId = -1;
        hipGetDevice(&deviceId);
        hipDeviceProp_t deviceProp;
        cutilSafeCall(hipGetDeviceProperties(&deviceProp, deviceId));
        printf("\n Using AMD device: %d %s\n\n", deviceId, deviceProp.name);
    } 
    else 
    {
        printf("\n Using CPU\n");
    }
}

void print_last_cuda_error()
{
    printf("ANY CUDA ERRORS? %s\n", hipGetErrorString(hipGetLastError()));
}










// **************************************** REMOVED FUNCTIONS ****************************************

// TODO: Try to do this in cuda
// make cell index of the particle to be removed as int_max,
// sort all arrays
// make sizes changed instead of resizing
// Could do this inside device itself
//****************************************
void oppic_finalize_particle_move_cuda(oppic_set set)
{
    cutilSafeCall(hipMemcpy(set->particle_statuses, set->particle_statuses_d, set->size * sizeof(int), hipMemcpyDeviceToHost));

    set->particle_remove_count = 0;

    for (int j = 0; j < set->size; j++)
    {
        if (set->particle_statuses[j] == OPP_NEED_REMOVE)
        {
            (set->particle_remove_count)++; // Could use a device variable and use atomicAdds inside device code
        }
    }   

    if (OP_DEBUG) printf("\toppic_finalize_particle_move_cuda set [%s] with particle_remove_count [%d]\n", set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0)
    {
        cutilSafeCall(hipFree(set->particle_statuses_d));
        set->particle_statuses_d = NULL;
        return;
    }

    oppic_download_particle_set(set); // TODO : Find a better way

    oppic_finalize_particle_move_core(set);

    oppic_upload_particle_set(set, true /* realloc the device pointers */);

    cutilSafeCall(hipFree(set->particle_statuses_d));
    set->particle_statuses_d = NULL;

    if (OP_DEBUG) printf("\toppic_finalize_particle_move_cuda set [%s] with new size [%d]\n", set->name, set->size);
}

// **************************************** ***************** ****************************************