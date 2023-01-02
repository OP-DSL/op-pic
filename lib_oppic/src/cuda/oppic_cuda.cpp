
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

#include <oppic_cuda.h>

//****************************************
void oppic_init(int argc, char **argv, int diags)
{
    oppic_init_core(argc, argv, diags);
    cutilDeviceInit(argc, argv);

    cutilSafeCall(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
    cutilSafeCall(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
}

//****************************************
void oppic_exit()
{
    oppic_cuda_exit();
    oppic_exit_core();
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

    op_cpHostToDevice((void **)&(map->map_d), (void **)&(temp_map), map->dim * set_size * sizeof(int));
    free(temp_map);

    return map;
}

//****************************************
oppic_dat oppic_decl_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name)
{
    return oppic_decl_dat_core(set, dim, type, size, data, name);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, int dim, const char *typ, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, idx, map, dim, typ, acc, map_with_cell_index);
}

//****************************************
oppic_arg oppic_arg_dat(oppic_dat dat, int idx, oppic_map map, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, idx, map, acc, map_with_cell_index);
}
oppic_arg oppic_arg_dat(oppic_dat dat, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(dat, acc, map_with_cell_index);
}
oppic_arg oppic_arg_dat(oppic_map map, oppic_access acc, bool map_with_cell_index)
{
    return oppic_arg_dat_core(map, acc, map_with_cell_index);
}

//****************************************
// template <class T> oppic_arg oppic_arg_gbl(T *data, int dim, char const *typ, oppic_access acc);
oppic_arg oppic_arg_gbl(double *data, int dim, char const *typ, oppic_access acc)
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
    return oppic_decl_particle_set_core(size, name, cells_set);
}
oppic_set oppic_decl_particle_set(char const *name, oppic_set cells_set)
{
    return oppic_decl_particle_set_core(name, cells_set);
}

//****************************************
oppic_dat oppic_decl_particle_dat(oppic_set set, int dim, char const *type, int size, char *data, char const *name, bool cell_index)
{
    return oppic_decl_particle_dat_core(set, dim, type, size, data, name, cell_index);
}

//****************************************
void oppic_increase_particle_count(oppic_set particles_set, const int num_particles_to_insert)
{
    oppic_increase_particle_count_core(particles_set, num_particles_to_insert);
}

//****************************************
void oppic_reset_num_particles_to_insert(oppic_set set)
{
    oppic_reset_num_particles_to_insert_core(set);
}

//****************************************
void oppic_mark_particle_to_remove(oppic_set set, int particle_index)
{
    oppic_mark_particle_to_remove_core(set, particle_index);
}

//****************************************
void oppic_remove_marked_particles_from_set(oppic_set set)
{
    oppic_remove_marked_particles_from_set_core(set);
}
void oppic_remove_marked_particles_from_set(oppic_set set, std::vector<int>& idx_to_remove)
{
    oppic_remove_marked_particles_from_set_core(set, idx_to_remove);
}

//****************************************
void oppic_particle_sort(oppic_set set)
{ TRACE_ME;

    oppic_particle_sort_core(set);
}

//****************************************






void __cudaSafeCall(cudaError_t err, const char *file, const int line) 
{
  if (cudaSuccess != err) {
    fprintf(stderr, "%s(%i) : cutilSafeCall() Runtime API error : %s.\n", file, line, cudaGetErrorString(err));
    exit(-1);
  }
}

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line) 
{
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err) {
    fprintf(stderr, "%s(%i) : cutilCheckMsg() error : %s : %s.\n", file, line, errorMessage, cudaGetErrorString(err));
    exit(-1);
  }
}

void op_mvHostToDevice(void **map, int size) 
{
  if (!OP_hybrid_gpu || size == 0)
    return;
  void *tmp;
  cutilSafeCall(cudaMalloc(&tmp, size));
  cutilSafeCall(cudaMemcpy(tmp, *map, size, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaDeviceSynchronize());
  free(*map);
  *map = tmp;
}

void op_cpHostToDevice(void **data_d, void **data_h, int size) 
{
  if (!OP_hybrid_gpu)
    return;
  if (*data_d != NULL) cutilSafeCall(cudaFree(*data_d));
  cutilSafeCall(cudaMalloc(data_d, size));
  cutilSafeCall(cudaMemcpy(*data_d, *data_h, size, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaDeviceSynchronize());
}

void cutilDeviceInit(int argc, char **argv) 
{
  (void)argc;
  (void)argv;
  int deviceCount;
  cutilSafeCall(cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    printf("cutil error: no devices supporting CUDA\n");
    exit(-1);
  }

  // Test we have access to a device
  float *test;
  cudaError_t err = cudaMalloc((void **)&test, sizeof(float));
  if (err != cudaSuccess) {
    OP_hybrid_gpu = 0;
  } else {
    OP_hybrid_gpu = 1;
  }
  if (OP_hybrid_gpu) {
    cudaFree(test);

    cutilSafeCall(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    int deviceId = -1;
    cudaGetDevice(&deviceId);
    cudaDeviceProp_t deviceProp;
    cutilSafeCall(cudaGetDeviceProperties(&deviceProp, deviceId));
    printf("\n Using CUDA device: %d %s\n", deviceId, deviceProp.name);
  } else {
    printf("\n Using CPU\n");
  }
}

void op_upload_dat(oppic_dat dat) {
  if (!OP_hybrid_gpu)
    return;
  size_t set_size = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
  if (strstr(dat->type, ":soa") != NULL || (OP_auto_soa && dat->dim > 1)) {
    char *temp_data = (char *)malloc(dat->size * set_size * sizeof(char));
    int element_size = dat->size / dat->dim;
    for (int i = 0; i < dat->dim; i++) {
      for (int j = 0; j < set_size; j++) {
        for (int c = 0; c < element_size; c++) {
          temp_data[element_size * i * set_size + element_size * j + c] =
              dat->data[dat->size * j + element_size * i + c];
        }
      }
    }
    cutilSafeCall(cudaMemcpy(dat->data_d, temp_data, set_size * dat->size,
                             cudaMemcpyHostToDevice));
    free(temp_data);
  } else {
    cutilSafeCall(cudaMemcpy(dat->data_d, dat->data, set_size * dat->size,
                             cudaMemcpyHostToDevice));
  }
}

void op_download_dat(oppic_dat dat) {
  if (!OP_hybrid_gpu)
    return;
  size_t set_size = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
  if (strstr(dat->type, ":soa") != NULL || (OP_auto_soa && dat->dim > 1)) {
    char *temp_data = (char *)malloc(dat->size * set_size * sizeof(char));
    cutilSafeCall(cudaMemcpy(temp_data, dat->data_d, set_size * dat->size,
                             cudaMemcpyDeviceToHost));
    int element_size = dat->size / dat->dim;
    for (int i = 0; i < dat->dim; i++) {
      for (int j = 0; j < set_size; j++) {
        for (int c = 0; c < element_size; c++) {
          dat->data[dat->size * j + element_size * i + c] =
              temp_data[element_size * i * set_size + element_size * j + c];
        }
      }
    }
    free(temp_data);
  } else {
    cutilSafeCall(cudaMemcpy(dat->data, dat->data_d, set_size * dat->size,
                             cudaMemcpyDeviceToHost));
  }
}

int op_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args) {
  for (int n = 0; n < nargs; n++)
    if (args[n].opt && args[n].argtype == OP_ARG_DAT &&
        args[n].dat->dirty_hd == Dirty::Host) {
      op_download_dat(args[n].dat);
      args[n].dat->dirty_hd = Dirty::NotDirty;
    }
  return set->size;
}

int op_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device){
  (void)device;
  return device == 1 ? op_mpi_halo_exchanges(set, nargs, args) : op_mpi_halo_exchanges_cuda(set, nargs, args);
}

void op_mpi_set_dirtybit(int nargs, oppic_arg *args) {
  for (int n = 0; n < nargs; n++) {
    if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
        (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
         args[n].acc == OP_RW)) {
      args[n].dat->dirty_hd = Dirty::Device;
    }
  }
}

int op_mpi_halo_exchanges_cuda(oppic_set set, int nargs, oppic_arg *args) { 
  for (int n = 0; n < nargs; n++)
  { 
    if (args[n].opt && args[n].argtype == OP_ARG_DAT &&
        args[n].dat->dirty_hd == Dirty::Device) { 
      op_upload_dat(args[n].dat);
      args[n].dat->dirty_hd = Dirty::NotDirty;
    }
  }
  return set->size;
}

void op_mpi_set_dirtybit_cuda(int nargs, oppic_arg *args) {
  for (int n = 0; n < nargs; n++) {
    if ((args[n].opt == 1) && (args[n].argtype == OP_ARG_DAT) &&
        (args[n].acc == OP_INC || args[n].acc == OP_WRITE ||
         args[n].acc == OP_RW)) {
      args[n].dat->dirty_hd = Dirty::Host;
    }
  }
}




void oppic_cuda_exit() {
  if (!OP_hybrid_gpu)
    return;

  for (auto& a : oppic_maps) {
    cutilSafeCall(cudaFree(a->map_d));
  }

  for (auto& a : oppic_dats) {
    cutilSafeCall(cudaFree(a->data_d));
  }
}

//****************************************
