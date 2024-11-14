

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

#include "opp_hip.h"
#ifdef USE_MPI
    #include "opp_hip_kernels.cpp"
#endif

enum UpDownType {
    None = 0,
    All = 1,
    AllHalo,
    OnlyExecute,
};

void opp_mv_halo_list_device();
void __opp_mpi_device_halo_exchange(opp_arg *arg, int exec_flag);
void __opp_mpi_device_halo_wait_all(opp_arg *arg);
void __opp_mpi_device_halo_wait_all(int nargs, opp_arg *args);

//*******************************************************************************
void opp_device_malloc(void **ptr, size_t size)
{
    if (*ptr != NULL) 
        cutilSafeCall(hipFree(*ptr));
    
    cutilSafeCall(hipMalloc(ptr, size));
}

//*******************************************************************************
void opp_device_free(void* ptr)
{
    cutilSafeCall(hipFree(ptr));
}

/*******************************************************************************
 * Main MPI halo creation routine
 *******************************************************************************/
void opp_halo_create() 
{
#ifdef USE_MPI
    __opp_halo_create();
    
    // The device arrays are dirty at this point, but opp_partition() routine calling this
    // makes all device arrays clean

    for (auto& set : opp_sets) 
    {
        if (set->is_particle) continue;

        for (auto& dat : opp_dats) 
        {
            if (dat->set->index == set->index)
            {
                if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat->dim > 1)) 
                {
                    opp_device_malloc((void **)&(dat->buffer_d_r),
                                (size_t)dat->size * (OPP_import_exec_list[set->index]->size +
                                                OPP_import_nonexec_list[set->index]->size));
                    if (OPP_DBG)
                        opp_printf("opp_halo_create", "buffer_d_r Alloc %zu bytes for %s dat", 
                            (size_t)dat->size * (OPP_import_exec_list[set->index]->size +
                            OPP_import_nonexec_list[set->index]->size), dat->name);
                } 

                opp_device_malloc((void **)&(dat->buffer_d),
                                (size_t)dat->size * (OPP_export_exec_list[set->index]->size +
                                            OPP_export_nonexec_list[set->index]->size)); 
                                            // + set_import_buffer_size[set->index]));
                
                if (OPP_DBG)
                    opp_printf("opp_halo_create", "buffer_d Alloc %zu bytes for %s dat", 
                        (size_t)dat->size * (OPP_export_exec_list[set->index]->size +
                        OPP_export_nonexec_list[set->index]->size), dat->name);
            }
        }
    }

    opp_mv_halo_list_device();
#endif
}

/*******************************************************************************
 * Routine to Clean-up all MPI halos(called at the end of an OP2 MPI application)
*******************************************************************************/
void opp_halo_destroy() 
{
#ifdef USE_MPI
    __opp_halo_destroy();

    if (OPP_hybrid_gpu) 
    {
        for (auto& dat : opp_dats) 
        {
            if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat->dim > 1)) 
            {
                cutilSafeCall(hipFree(dat->buffer_d_r));
            }
            cutilSafeCall(hipFree(dat->buffer_d));
        }

        for (int i = 0; i < (int)opp_sets.size(); i++) 
        {
            if (export_exec_list_d[i] != NULL)
                cutilSafeCall(hipFree(export_exec_list_d[i]));
            if (export_nonexec_list_d[i] != NULL)
                cutilSafeCall(hipFree(export_nonexec_list_d[i]));
        }
    }
#endif
}

/*******************************************************************************
 * Routine to exchange MPI halos of the all the args
*******************************************************************************/
int opp_mpi_halo_exchanges(opp_set set, int nargs, opp_arg *args) 
{

    if (OPP_DBG) opp_printf("opp_mpi_halo_exchanges", "START");

    int size = set->size;

#ifdef USE_MPI
    bool direct_flag = true;

    // check if this is a direct loop
    for (int n = 0; n < nargs; n++)
        if (args[n].opt && args[n].argtype == OPP_ARG_DAT && args[n].idx != -1)
            direct_flag = false;

    // return set size if it is a direct loop
    if (direct_flag)
    {
        if (OPP_DBG) opp_printf("opp_mpi_halo_exchanges", "This is a direct loop");
        return size;
    }   

    if (OPP_DBG) 
        opp_printf("opp_mpi_halo_exchanges", "This is a in-direct loop");

    int exec_flag = 0;
    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].idx != -1 && args[n].acc != OPP_READ) 
        {
            size = set->size + set->exec_size;
            exec_flag = 1;
        }
    }

    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].argtype == OPP_ARG_DAT && (!args[n].dat->set->is_particle) && 
            (args[n].dat->dirtybit == 1)) 
        {
            bool already_done = false;

            // Check if dat was already done within these args
            for (int m = 0; m < n; m++) 
            {
                if (args[n].dat == args[m].dat)
                    already_done = true;
            }

            if (!already_done)
            {
                if (OPP_DBG) 
                    opp_printf("opp_mpi_halo_exchanges OLD", "opp_exchange_halo for dat [%s] exec_flag %d", 
                        args[n].dat->name, exec_flag);
                __opp_mpi_host_halo_exchange(&args[n], exec_flag);
            }
        }
    }  
#endif

    return size;
}

// if multiple of halo exchanges are done for different kernels at once, this may get currupted
DeviceType opp_current_device = Device_CPU;

struct opp_HaloExInfo {
    bool skip = false;
    UpDownType download = UpDownType::None;
    bool HaloEx = false;
    UpDownType upload = UpDownType::None;
};

std::vector<opp_HaloExInfo> haloExInfo;
int current_device = Device_GPU;

/*
RUN_ON	LOOPP_TYPE	DirtyBit	DirtyHD	    Download	HaloEx	Upload	| Set_DirtyBit	Set_DirtyHD
----------------------------------------------------------------------------------------------------
DEVICE	DIRECT LOOP	    0	    Not Dirty		0	    0	    0	    |   0	        Not Dirty
		                1	    Not Dirty		0	    0	    0	    |   1	        Not Dirty
		                0	    DeviceDirty		0	    0	    1	    |   0	        Not Dirty
		                1	    DeviceDirty		0	    0	    1	    |   1	        Not Dirty
		                0	    Host Dirty		0	    0	    0	    |   0	        Host Dirty
		                1	    Host Dirty		0	    0	    0	    |   1	        Host Dirty
       
	    INDIRECT LOOP	0	    Not Dirty		0	    0	    0	    |   0	        Not Dirty
		                1	    Not Dirty		0	    1	    1	    |   0	        Not Dirty
		                0	    DeviceDirty		0	    0	    1	    |   0	        Not Dirty
		                1	    DeviceDirty		0	    1	    1	    |   0	        Not Dirty
		                0	    Host Dirty		0	    0	    0	    |   0	        Host Dirty
		                1	    Host Dirty		1	    1	    1	    |   0	        Not Dirty
       
HOST	DIRECT LOOP	    0	    Not Dirty		0	    0	    0	    |   0	        Not Dirty
                        1	    Not Dirty		0	    0	    0	    |   1	        Not Dirty
                        0	    DeviceDirty		0	    0	    0	    |   0	        DeviceDirty
                        1	    DeviceDirty		0	    0	    0	    |   1	        DeviceDirty
                        0	    Host Dirty		1	    0	    0	    |   0	        Not Dirty
                        1	    Host Dirty		1	    0	    0	    |   1	        Not Dirty
                           
        INDIRECT LOOP	0	    Not Dirty		0	    0	    0	    |   0	        Not Dirty
                        1	    Not Dirty		0	    1	    0	    |   0	        Not Dirty
                        0	    DeviceDirty		0	    0	    0	    |   0	        DeviceDirty
                        1	    DeviceDirty		0	    1	    0	    |   0	        DeviceDirty
                        0	    Host Dirty		1	    0	    0	    |   0	        Not Dirty
                        1	    Host Dirty		1	    1	    0	    |   0	        Not Dirty
*/

void markExInfo(opp_arg& arg, DeviceType device, bool direct_loop, opp_HaloExInfo& exInfo) {

    // return if the arg is not for a dat or it is not read somehow
    if (!arg.opt || arg.argtype != OPP_ARG_DAT || !(arg.acc == OPP_READ || arg.acc == OPP_RW))
        return;
    
    // set halo ex flag, only if loop is indirect and DirtyBit is set
    if (!arg.dat->set->is_particle && !direct_loop && arg.dat->dirtybit == 1) {
        exInfo.HaloEx = true;
    }

    // set download flag, only if a dat is to be read somehow and 
    // host is dirty while operating on CPU or host is dirty while operating on GPU but need halo ex
    if (device == Device_CPU && arg.dat->dirty_hd == Dirty::Host) {
        exInfo.download = UpDownType::All;
    }
    // else if (device == Device_GPU && !direct_loop && arg.dat->dirty_hd == Dirty::Host 
    //             && arg.dat->dirtybit == 1 && OPP_comm_size > 1))
    //     exInfo.download = UpDownType::AllHalo;

    // set upload flag, only if operating on GPU and if device is dirty and/or need a halo Ex
    if (device == Device_GPU)
    {
        if (arg.dat->dirty_hd == Dirty::Device)
            exInfo.upload = UpDownType::All;
        // else if (!direct_loop && arg.dat->dirtybit == 1  && OPP_comm_size > 1)
        //     exInfo.upload = UpDownType::AllHalo;
    } 

}

void changeDatFlags(opp_arg& arg, DeviceType device, bool direct_loop) {

    if (!arg.opt || arg.argtype != OPP_ARG_DAT || !(arg.acc == OPP_READ || arg.acc == OPP_RW))
        return;

    if (device == Device_GPU && arg.dat->dirty_hd == Dirty::Host) // running on device so host wont get fully updated
            arg.dat->dirty_hd = Dirty::Host;
    else if (device == Device_CPU && arg.dat->dirty_hd == Dirty::Device) // running on host so device wont get updated
            arg.dat->dirty_hd = Dirty::Device;
    else
            arg.dat->dirty_hd = Dirty::NotDirty;

    // if not a direct loop, then halo ex will happen, so mark dirtybit to 0
    if (!direct_loop)
        arg.dat->dirtybit = 0;
}

void generateHaloExchangeInfo(opp_set iter_set, int nargs, opp_arg *args, DeviceType device) {

    haloExInfo.clear();
    bool direct_loop = true;

    // check whether the loop is a direct loop or not
    for (int n = 0; n < nargs; n++)
        if (args[n].opt && args[n].argtype == OPP_ARG_DAT && args[n].idx != -1 && 
            ((!iter_set->is_particle && !args[n].dat->set->is_particle) || 
            (iter_set->is_particle && args[n].dat->set->index != iter_set->cells_set->index)))
                direct_loop = false;

    for (int n = 0; n < nargs; n++) {

        opp_HaloExInfo exInfo;

        bool already_done = false;
        for (int m = 0; m < n; m++) { // Check if dat was already done within these args
            if (args[n].dat == args[m].dat)
                already_done = true;
        }

        if (!already_done) {
            markExInfo(args[n], device, direct_loop, exInfo);
        }
        else {
            exInfo.skip = true;
        }

        haloExInfo.push_back(exInfo);
    }

    for (int n = 0; n < nargs; n++) {

        if (!haloExInfo[n].skip) {
            changeDatFlags(args[n], device, direct_loop);
        }
    }
}


int opp_mpi_halo_exchanges_grouped(opp_set set, int nargs, opp_arg *args, DeviceType device)
{
    current_device = device;
    int size = set->size;
    
    generateHaloExchangeInfo(set, nargs, args, device);

    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].download == UpDownType::All) {
            if (OPP_DBG) opp_printf("halo_ex", "downloading %s", args[n].dat->name);
            opp_download_dat(args[n].dat);
        }
    }

#ifdef USE_MPI
    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].HaloEx) {
            if (device == Device_CPU)
            {
                if (OPP_DBG) opp_printf("halo_ex", "__opp_mpi_host_halo_exchange for dat [%s] exec_flag %d", args[n].dat->name, 1);
                __opp_mpi_host_halo_exchange(&args[n], 1);   
            }
            else // Device_GPU
            {
                if (OPP_DBG) opp_printf("halo_ex", "__opp_mpi_device_halo_exchange for dat [%s] exec_flag %d", args[n].dat->name, 1);
                __opp_mpi_device_halo_exchange(&args[n], 1);
            }
            size = set->size + set->exec_size;
        }
    }
#endif

    return size;
}

void opp_mpi_force_halo_update_if_dirty(opp_set set, std::vector<opp_dat> dats, DeviceType device) {
    
    const int nargs = (int)dats.size();
    std::vector<opp_arg> args(nargs);

    for (int i = 0; i < nargs; i++) {
        args[i] = opp_arg_dat(dats[i], OPP_READ);
        args[i].idx = 2; // HACK to forcefully make halos to download
    }

    opp_mpi_halo_exchanges_grouped(set, nargs, args.data(), device);
    opp_mpi_halo_wait_all(nargs, args.data());
}

void opp_mpi_halo_wait_all(int nargs, opp_arg *args)
{
#ifdef USE_MPI
    // Finish the halo exchange
    if (current_device == Device_CPU)
    {
        __opp_mpi_host_halo_wait_all(nargs, args); // Halo exchange only for mesh dats
    }
    else // Device_GPU
    {
        __opp_mpi_device_halo_wait_all(nargs, args); // Halo exchange only for mesh dats
    }
#endif

    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].upload == UpDownType::All) {
            if (OPP_DBG) opp_printf("halo_ex", "uploading %s", args[n].dat->name);
            opp_upload_dat(args[n].dat); // TODO : ideally need to upload only the halo regions
        }
    }
}

void __opp_mpi_device_halo_exchange(opp_arg *arg, int exec_flag) 
{
#ifdef USE_MPI
    opp_dat dat = arg->dat;

    if (arg->sent == 1) 
    {
        // opp_printf("opp_mpi_halo_exchange_dev", "Error: Halo exchange already in flight for dat %s", dat->name);
        fflush(stdout);
        opp_abort("opp_mpi_halo_exchange_dev");
    }

    arg->sent = 0; // reset flag
    // need to exchange both direct and indirect data sets if they are dirty

    if (OPP_DBG) 
        opp_printf("opp_mpi_halo_exchange_dev", "Exchanging Halo of data array [%s]", dat->name);

    halo_list imp_exec_list = OPP_import_exec_list[dat->set->index];
    halo_list imp_nonexec_list = OPP_import_nonexec_list[dat->set->index];

    halo_list exp_exec_list = OPP_export_exec_list[dat->set->index];
    halo_list exp_nonexec_list = OPP_export_nonexec_list[dat->set->index];

    //-------first exchange exec elements related to this data array--------

    // sanity checks
    if (compare_sets(imp_exec_list->set, dat->set) == 0) {
        printf("Error: Import list and set mismatch\n");
        MPI_Abort(OPP_MPI_WORLD, 2);
    }
    if (compare_sets(exp_exec_list->set, dat->set) == 0) {
        printf("Error: Export list and set mismatch\n");
        MPI_Abort(OPP_MPI_WORLD, 2);
    }

    gather_data_to_buffer(*arg, exp_exec_list, exp_nonexec_list);

    char *outptr_exec = NULL;
    char *outptr_nonexec = NULL;
    if (OPP_gpu_direct) {
        outptr_exec = arg->dat->buffer_d;
        outptr_nonexec =
            arg->dat->buffer_d + exp_exec_list->size * arg->dat->size;
        cutilSafeCall(hipDeviceSynchronize());
    } else {
        cutilSafeCall(hipMemcpy(
            ((op_mpi_buffer)(dat->mpi_buffer))->buf_exec, arg->dat->buffer_d,
            exp_exec_list->size * arg->dat->size, hipMemcpyDeviceToHost));

        cutilSafeCall(hipMemcpy(
            ((op_mpi_buffer)(dat->mpi_buffer))->buf_nonexec,
            arg->dat->buffer_d + exp_exec_list->size * arg->dat->size,
            exp_nonexec_list->size * arg->dat->size, hipMemcpyDeviceToHost));

        cutilSafeCall(hipDeviceSynchronize());
        outptr_exec = ((op_mpi_buffer)(dat->mpi_buffer))->buf_exec;
        outptr_nonexec = ((op_mpi_buffer)(dat->mpi_buffer))->buf_nonexec;
    }

    for (int i = 0; i < exp_exec_list->ranks_size; i++) {
        MPI_Isend(&outptr_exec[exp_exec_list->disps[i] * dat->size],
                    dat->size * exp_exec_list->sizes[i], MPI_CHAR,
                    exp_exec_list->ranks[i], dat->index, OPP_MPI_WORLD,
                    &((op_mpi_buffer)(dat->mpi_buffer))
                        ->s_req[((op_mpi_buffer)(dat->mpi_buffer))->s_num_req++]);
    }

    int init = dat->set->size * dat->size;
    char *ptr = NULL;
    for (int i = 0; i < imp_exec_list->ranks_size; i++) {
        ptr = OPP_gpu_direct
                    ? &(dat->data_d[init + imp_exec_list->disps[i] * dat->size])
                    : &(dat->data[init + imp_exec_list->disps[i] * dat->size]);
        if (OPP_gpu_direct && (strstr(arg->dat->type, ":soa") != NULL ||
                                (OPP_auto_soa && arg->dat->dim > 1)))
            ptr = dat->buffer_d_r + imp_exec_list->disps[i] * dat->size;
        MPI_Irecv(ptr, dat->size * imp_exec_list->sizes[i], MPI_CHAR,
                    imp_exec_list->ranks[i], dat->index, OPP_MPI_WORLD,
                    &((op_mpi_buffer)(dat->mpi_buffer))
                        ->r_req[((op_mpi_buffer)(dat->mpi_buffer))->r_num_req++]);
    }

    //-----second exchange nonexec elements related to this data array------
    // sanity checks
    if (compare_sets(imp_nonexec_list->set, dat->set) == 0) {
        printf("Error: Non-Import list and set mismatch");
        MPI_Abort(OPP_MPI_WORLD, 2);
    }
    if (compare_sets(exp_nonexec_list->set, dat->set) == 0) {
        printf("Error: Non-Export list and set mismatch");
        MPI_Abort(OPP_MPI_WORLD, 2);
    }

    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) {
        MPI_Isend(&outptr_nonexec[exp_nonexec_list->disps[i] * dat->size],
                    dat->size * exp_nonexec_list->sizes[i], MPI_CHAR,
                    exp_nonexec_list->ranks[i], dat->index, OPP_MPI_WORLD,
                    &((op_mpi_buffer)(dat->mpi_buffer))
                        ->s_req[((op_mpi_buffer)(dat->mpi_buffer))->s_num_req++]);
    }

    int nonexec_init = (dat->set->size + imp_exec_list->size) * dat->size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) {
        ptr = OPP_gpu_direct
                    ? &(dat->data_d[nonexec_init +
                                    imp_nonexec_list->disps[i] * dat->size])
                    : &(dat->data[nonexec_init +
                                imp_nonexec_list->disps[i] * dat->size]);
        if (OPP_gpu_direct && (strstr(arg->dat->type, ":soa") != NULL ||
                                (OPP_auto_soa && arg->dat->dim > 1)))
            ptr = dat->buffer_d_r +
                (imp_exec_list->size + imp_exec_list->disps[i]) * dat->size;
        MPI_Irecv(ptr, dat->size * imp_nonexec_list->sizes[i], MPI_CHAR,
                    imp_nonexec_list->ranks[i], dat->index, OPP_MPI_WORLD,
                    &((op_mpi_buffer)(dat->mpi_buffer))
                        ->r_req[((op_mpi_buffer)(dat->mpi_buffer))->r_num_req++]);
    }

    // clear dirty bit
    dat->dirtybit = 0;
    arg->sent = 1;
#endif
}



void __opp_mpi_device_halo_wait_all(opp_arg *arg) 
{
#ifdef USE_MPI
    opp_dat dat = arg->dat;

    MPI_Waitall(((op_mpi_buffer)(dat->mpi_buffer))->s_num_req,
                ((op_mpi_buffer)(dat->mpi_buffer))->s_req, MPI_STATUSES_IGNORE);
    MPI_Waitall(((op_mpi_buffer)(dat->mpi_buffer))->r_num_req,
                ((op_mpi_buffer)(dat->mpi_buffer))->r_req, MPI_STATUSES_IGNORE);

    ((op_mpi_buffer)(dat->mpi_buffer))->s_num_req = 0;
    ((op_mpi_buffer)(dat->mpi_buffer))->r_num_req = 0;

    if (OPP_gpu_direct == 0) 
    {

        if (strstr(arg->dat->type, ":soa") != NULL || (OPP_auto_soa && arg->dat->dim > 1)) 
        {
            int init = dat->set->size * dat->size;
            int size = (dat->set->exec_size + dat->set->nonexec_size) * dat->size;
            cutilSafeCall(hipMemcpyAsync(dat->buffer_d_r, dat->data + init, size,
                                            hipMemcpyHostToDevice, 0));
            scatter_data_from_buffer(*arg);
        } 
        else 
        {
            int init = dat->set->size * dat->size;
            cutilSafeCall(hipMemcpyAsync(dat->data_d + init, dat->data + init,
                            (OPP_import_exec_list[dat->set->index]->size +
                            OPP_import_nonexec_list[dat->set->index]->size) *
                                arg->dat->size,
                            hipMemcpyHostToDevice, 0));
        }
    } 
    else if (strstr(arg->dat->type, ":soa") != NULL || (OPP_auto_soa && arg->dat->dim > 1))
    {
        scatter_data_from_buffer(*arg);
    }

    arg->sent = 2; // set flag to indicate completed comm
#endif
}


void __opp_mpi_device_halo_wait_all(int nargs, opp_arg *args) 
{
    if (OPP_DBG) opp_printf("__opp_mpi_device_halo_wait_all", "START");

    opp_profiler->startMpiComm("", opp::OPP_Mesh);

    for (int n = 0; n < nargs; n++) 
    {
        opp_arg *arg = &args[n];

        if (arg->opt && arg->argtype == OPP_ARG_DAT && arg->sent == 1) 
        {
            __opp_mpi_device_halo_wait_all(arg);
        }
    }

    opp_profiler->endMpiComm("", opp::OPP_Mesh);

    if (OPP_DBG) opp_printf("__opp_mpi_device_halo_wait_all", "END");
}





void opp_mv_halo_list_device() 
{
    if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "START");

#ifdef USE_MPI

    if (export_exec_list_d != NULL) 
    {
        for (int s = 0; s < (int)opp_sets.size(); s++)
            if (export_exec_list_d[opp_sets[s]->index] != NULL)
                cutilSafeCall(hipFree(export_exec_list_d[opp_sets[s]->index]));
        free(export_exec_list_d);
    }
    export_exec_list_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    {
        opp_set set = opp_sets[s];
        export_exec_list_d[set->index] = NULL;

        if (set->is_particle) continue;

        opp_cpHostToDevice((void **)&(export_exec_list_d[set->index]),
                        (void **)&(OPP_export_exec_list[set->index]->list),
                        OPP_export_exec_list[set->index]->size * sizeof(int),
                        OPP_export_exec_list[set->index]->size * sizeof(int), true);

        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "export_exec_list_d Alloc %zu bytes for set %s", 
            OPP_export_exec_list[set->index]->size * sizeof(int), set->name);
    }

    if (export_nonexec_list_d != NULL) 
    {
        for (int s = 0; s < (int)opp_sets.size(); s++)
            if (export_nonexec_list_d[opp_sets[s]->index] != NULL)
                cutilSafeCall(hipFree(export_nonexec_list_d[opp_sets[s]->index]));
        free(export_nonexec_list_d);
    }
    export_nonexec_list_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    {
        opp_set set = opp_sets[s];
        export_nonexec_list_d[set->index] = NULL;
        
        if (set->is_particle) continue;

        opp_cpHostToDevice((void **)&(export_nonexec_list_d[set->index]),
                        (void **)&(OPP_export_nonexec_list[set->index]->list),
                        OPP_export_nonexec_list[set->index]->size * sizeof(int),
                        OPP_export_nonexec_list[set->index]->size * sizeof(int), true);
        
        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "export_nonexec_list_d Alloc %zu bytes for set %s", 
            OPP_export_nonexec_list[set->index]->size * sizeof(int), set->name);
    }

    //for grouped, we need the disps array on device too
    if (export_exec_list_disps_d != NULL) 
    {
        for (int s = 0; s < (int)opp_sets.size(); s++)
            if (export_exec_list_disps_d[opp_sets[s]->index] != NULL)
                cutilSafeCall(hipFree(export_exec_list_disps_d[opp_sets[s]->index]));
        free(export_exec_list_disps_d);
    }
    export_exec_list_disps_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];
        export_exec_list_disps_d[set->index] = NULL;

        if (set->is_particle) continue;

        //make sure end size is there too
        OPP_export_exec_list[set->index]
            ->disps[OPP_export_exec_list[set->index]->ranks_size] =
            OPP_export_exec_list[set->index]->ranks_size == 0
                ? 0
                : OPP_export_exec_list[set->index]
                        ->disps[OPP_export_exec_list[set->index]->ranks_size - 1] +
                    OPP_export_exec_list[set->index]
                        ->sizes[OPP_export_exec_list[set->index]->ranks_size - 1];
        opp_cpHostToDevice((void **)&(export_exec_list_disps_d[set->index]),
                        (void **)&(OPP_export_exec_list[set->index]->disps),
                        (OPP_export_exec_list[set->index]->ranks_size+1) * sizeof(int),
                        (OPP_export_exec_list[set->index]->ranks_size+1) * sizeof(int), true);
        
        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "export_exec_list_disps_d Alloc %zu bytes for set %s", 
            (OPP_export_exec_list[set->index]->ranks_size+1) * sizeof(int), set->name);
    }

    if (export_nonexec_list_disps_d != NULL) 
    {
        for (int s = 0; s < (int)opp_sets.size(); s++)
            if (export_nonexec_list_disps_d[opp_sets[s]->index] != NULL)
                cutilSafeCall(hipFree(export_nonexec_list_disps_d[opp_sets[s]->index]));
        free(export_nonexec_list_disps_d);
    }
    export_nonexec_list_disps_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++) { // for each set
        opp_set set = opp_sets[s];
        export_nonexec_list_disps_d[set->index] = NULL;

        if (set->is_particle) continue;

        //make sure end size is there too
        OPP_export_nonexec_list[set->index]
            ->disps[OPP_export_nonexec_list[set->index]->ranks_size] =
            OPP_export_nonexec_list[set->index]->ranks_size == 0
                ? 0
                : OPP_export_nonexec_list[set->index]
                        ->disps[OPP_export_nonexec_list[set->index]->ranks_size -
                                1] +
                    OPP_export_nonexec_list[set->index]
                        ->sizes[OPP_export_nonexec_list[set->index]->ranks_size -
                                1];
        opp_cpHostToDevice((void **)&(export_nonexec_list_disps_d[set->index]),
                        (void **)&(OPP_export_nonexec_list[set->index]->disps),
                        (OPP_export_nonexec_list[set->index]->ranks_size+1) * sizeof(int),
                        (OPP_export_nonexec_list[set->index]->ranks_size+1) * sizeof(int), true);
        
        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "export_nonexec_list_disps_d Alloc %zu bytes for set %s", 
            (OPP_export_nonexec_list[set->index]->ranks_size+1) * sizeof(int), set->name);
    }

    if (import_exec_list_disps_d != NULL) {
        for (int s = 0; s < (int)opp_sets.size(); s++)
        if (import_exec_list_disps_d[opp_sets[s]->index] != NULL)
            cutilSafeCall(hipFree(import_exec_list_disps_d[opp_sets[s]->index]));
        free(import_exec_list_disps_d);
    }
    import_exec_list_disps_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];
        import_exec_list_disps_d[set->index] = NULL;

        if (set->is_particle) continue;

        //make sure end size is there too
        OPP_import_exec_list[set->index]
            ->disps[OPP_import_exec_list[set->index]->ranks_size] =
            OPP_import_exec_list[set->index]->ranks_size == 0
                ? 0
                : OPP_import_exec_list[set->index]
                        ->disps[OPP_import_exec_list[set->index]->ranks_size - 1] +
                    OPP_import_exec_list[set->index]
                        ->sizes[OPP_import_exec_list[set->index]->ranks_size - 1];
        opp_cpHostToDevice((void **)&(import_exec_list_disps_d[set->index]),
                        (void **)&(OPP_import_exec_list[set->index]->disps),
                        (OPP_import_exec_list[set->index]->ranks_size+1) * sizeof(int),
                        (OPP_import_exec_list[set->index]->ranks_size+1) * sizeof(int), true);

        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "import_exec_list_disps_d Alloc %zu bytes for set %s", 
            (OPP_import_exec_list[set->index]->ranks_size+1) * sizeof(int), set->name);
    }

    if (import_nonexec_list_disps_d != NULL) 
    {
        for (int s = 0; s < (int)opp_sets.size(); s++)
            if (import_nonexec_list_disps_d[opp_sets[s]->index] != NULL)
                cutilSafeCall(hipFree(import_nonexec_list_disps_d[opp_sets[s]->index]));
        free(import_nonexec_list_disps_d);
    }
    import_nonexec_list_disps_d = (int **)malloc(sizeof(int *) * (int)opp_sets.size());

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];
        import_nonexec_list_disps_d[set->index] = NULL;

        if (set->is_particle) continue;
        
        //make sure end size is there too
        OPP_import_nonexec_list[set->index]
            ->disps[OPP_import_nonexec_list[set->index]->ranks_size] =
            OPP_import_nonexec_list[set->index]->ranks_size == 0
                ? 0
                : OPP_import_nonexec_list[set->index]
                        ->disps[OPP_import_nonexec_list[set->index]->ranks_size -
                                1] +
                    OPP_import_nonexec_list[set->index]
                        ->sizes[OPP_import_nonexec_list[set->index]->ranks_size -
                                1];
        opp_cpHostToDevice((void **)&(import_nonexec_list_disps_d[set->index]),
                        (void **)&(OPP_import_nonexec_list[set->index]->disps),
                        (OPP_import_nonexec_list[set->index]->ranks_size+1) * sizeof(int),
                        (OPP_import_nonexec_list[set->index]->ranks_size+1) * sizeof(int), true);

        if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "import_nonexec_list_disps_d Alloc %zu bytes for set %s", 
            (OPP_import_nonexec_list[set->index]->ranks_size+1) * sizeof(int), set->name);
    }

    if (OPP_DBG) opp_printf("opp_mv_halo_list_device", "END");
#endif
}