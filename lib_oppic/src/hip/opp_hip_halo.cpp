

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

/*******************************************************************************
 * Main MPI halo creation routine
 *******************************************************************************/
void opp_halo_create() 
{
#ifdef USE_MPI
    __opp_halo_create();
    
    // The device arrays are dirty at this point, but opp_partition() routine calling this
    // makes all device arrays clean
#endif
}

/*******************************************************************************
 * Routine to Clean-up all MPI halos(called at the end of an OP2 MPI application)
*******************************************************************************/
void opp_halo_destroy() 
{
#ifdef USE_MPI
    __opp_halo_destroy();
#endif
}

/*******************************************************************************
 * Routine to exchange MPI halos of the all the args
*******************************************************************************/
int opp_mpi_halo_exchanges(oppic_set set, int nargs, oppic_arg *args) 
{

    if (OP_DEBUG) opp_printf("opp_mpi_halo_exchanges", "START");

    int size = set->size;

#ifdef USE_MPI
    bool direct_flag = true;

    // check if this is a direct loop
    for (int n = 0; n < nargs; n++)
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].idx != -1)
            direct_flag = false;

    // return set size if it is a direct loop
    if (direct_flag)
    {
        if (OP_DEBUG) opp_printf("opp_mpi_halo_exchanges", "This is a direct loop");
        return size;
    }   

    if (OP_DEBUG) 
        opp_printf("opp_mpi_halo_exchanges", "This is a in-direct loop");

    int exec_flag = 0;
    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].idx != -1 && args[n].acc != OP_READ) 
        {
            size = set->size + set->exec_size;
            exec_flag = 1;
        }
    }

    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && (!args[n].dat->set->is_particle) && 
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
                if (OP_DEBUG) 
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
    bool download = false;
    bool HaloEx = false;
    bool upload = false;
};

std::vector<opp_HaloExInfo> haloExInfo;

/*
RUN_ON	LOOP_TYPE	DirtyBit	DirtyHD	    Download	HaloEx	Upload	| Set_DirtyBit	Set_DirtyHD
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

void markExInfo(oppic_arg& arg, DeviceType device, bool direct_loop, opp_HaloExInfo& exInfo) {

    // return if the arg is not for a dat or it is not read somehow
    if (!arg.opt || arg.argtype != OP_ARG_DAT || !(arg.acc == OP_READ || arg.acc == OP_RW))
        return;
    
    // set halo ex flag, only if loop is indirect and DirtyBit is set
    if (!arg.dat->set->is_particle && !direct_loop && arg.dat->dirtybit == 1) {
        exInfo.HaloEx = true;
    }

    // set download flag, only if a dat is to be read somehow and 
    // host is dirty while operating on CPU or host is dirty while operating on GPU but need halo ex
    if ((device == Device_CPU && arg.dat->dirty_hd == Dirty::Host) ||
        (device == Device_GPU && !direct_loop && arg.dat->dirty_hd == Dirty::Host && arg.dat->dirtybit == 1))
            exInfo.download = true;

    // set upload flag, only if operating on GPU and if device is dirty and/or need a halo Ex
    if (device == Device_GPU && (arg.dat->dirty_hd == Dirty::Device || (!direct_loop && arg.dat->dirtybit == 1)))
        exInfo.upload = true;
}

void changeDatFlags(oppic_arg& arg, DeviceType device, bool direct_loop) {

    if (!arg.opt || arg.argtype != OP_ARG_DAT || !(arg.acc == OP_READ || arg.acc == OP_RW))
        return;

    if (device == Device_GPU && arg.dat->dirty_hd == Dirty::Host && // halo exchange wont occur in below condition
        (direct_loop || (!direct_loop && arg.dat->dirtybit == 0)))
            arg.dat->dirty_hd = Dirty::Host;
    else if (device == Device_CPU && arg.dat->dirty_hd == Dirty::Device) // running on host so device wont get updated
            arg.dat->dirty_hd = Dirty::Device;
    else
            arg.dat->dirty_hd = Dirty::NotDirty;

    // if not a direct loop, then halo ex will happen, so mark dirtybit to 0
    if (!direct_loop)
        arg.dat->dirtybit = 0;
}

void generateHaloExchangeInfo(opp_set iter_set, int nargs, oppic_arg *args, DeviceType device) {

    haloExInfo.clear();
    bool direct_loop = true;

    // check whether the loop is a direct loop or not
    for (int n = 0; n < nargs; n++)
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].idx != -1 && 
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


int opp_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device)
{
    int size = set->size;
    
    generateHaloExchangeInfo(set, nargs, args, device);

    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].download) {
            // opp_printf("halo_ex", "downloading %s", args[n].dat->name);
            opp_download_dat(args[n].dat);
        }
    }

#ifdef USE_MPI
    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].HaloEx) {
            // opp_printf("halo_ex", "opp_exchange_halo for dat [%s] exec_flag %d", args[n].dat->name, 1);
            __opp_mpi_host_halo_exchange(&args[n], 1);
            size = set->size + set->exec_size;
        }
    }
#endif

    return size;
}

void opp_mpi_halo_wait_all(int nargs, oppic_arg *args)
{
#ifdef USE_MPI
    // Finish the halo exchange
    __opp_mpi_host_halo_wait_all(nargs, args); // Halo exchange only for mesh dats
#endif

    for (int n = 0; n < nargs; n++)
    {
        if (!haloExInfo[n].skip && haloExInfo[n].upload) {
            // opp_printf("halo_ex", "uploading %s", args[n].dat->name);
            opp_upload_dat(args[n].dat); // TODO : ideally need to upload only the halo regions
        }
    }
}


// /*******************************************************************************
//  * Routine to start exchanging halos for the args
// *******************************************************************************/
// int opp_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device)
// {
//     opp_current_device = device;

//     for (int n = 0; n < nargs; n++)
//     {
//         bool already_done = false;

//         // Check if dat reduction was already done within these args
//         for (int m = 0; m < n; m++) 
//         {
//             if (args[n].dat == args[m].dat)
//                 already_done = true;
//         }

//         if (!already_done && args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == Dirty::Host && 
//             (!args[n].dat->set->is_particle || opp_current_device == Device_CPU)) 
//         {
//             opp_printf("HALO OLD", "downloading %s", args[n].dat->name);
//             opp_download_dat(args[n].dat);
//             if (args[n].dat->dirty_hd == Dirty::Host)
//                 args[n].dat->dirty_hd = Dirty::NotDirty;
//         }
//     }

// #ifdef USE_MPI
//     return opp_mpi_halo_exchanges(set, nargs, args); // Halo exchange only for mesh dats
// #else
//     return set->size;
// #endif
// }

// /*******************************************************************************
//  * Routine to complete exchanging halos for the args
// *******************************************************************************/
// void opp_mpi_halo_wait_all(int nargs, oppic_arg *args)
// {
// #ifdef USE_MPI
//     // Finish the halo exchange
//     __opp_mpi_host_halo_wait_all(nargs, args); // Halo exchange only for mesh dats
// #endif

//     for (int n = 0; n < nargs; n++)
//     { 
//         bool already_done = false;

//         // Check if dat reduction was already done within these args
//         for (int m = 0; m < n; m++) 
//         {
//             if (args[n].dat == args[m].dat)
//                 already_done = true;
//         }

//         if (!already_done && args[n].opt && args[n].argtype == OP_ARG_DAT && opp_current_device == Device_GPU &&
//                 (args[n].dat->dirty_hd == Dirty::Device || !args[n].dat->set->is_particle)) 
//         { 
//             opp_printf("HALO OLD", "uploading %s", args[n].dat->name);
//             opp_upload_dat(args[n].dat); // TODO : ideally need to upload only the halo regions
//             if (args[n].dat->dirty_hd == Dirty::Device)
//                 args[n].dat->dirty_hd = Dirty::NotDirty;          
//         }
//     }        

//     for (int n = 0; n < nargs; n++) // This is wrong - change
//     { 
//         if (args[n].opt && args[n].argtype == OP_ARG_DAT && !args[n].dat->set->is_particle) 
//         {
//             if (opp_current_device == Device_CPU) 
//                 args[n].dat->dirty_hd = Dirty::Device;
//             else
//                 args[n].dat->dirty_hd = Dirty::NotDirty;
//         }
//     }
// }

