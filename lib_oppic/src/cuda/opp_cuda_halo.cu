

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

#include <opp_cuda.h>

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
                    opp_printf("opp_mpi_halo_exchanges", "opp_exchange_halo for dat [%s] exec_flag %d", 
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

/*******************************************************************************
 * Routine to start exchanging halos for the args
*******************************************************************************/
int opp_mpi_halo_exchanges_grouped(oppic_set set, int nargs, oppic_arg *args, DeviceType device)
{
    opp_current_device = device;

    for (int n = 0; n < nargs; n++)
    {
        bool already_done = false;

        // Check if dat reduction was already done within these args
        for (int m = 0; m < n; m++) 
        {
            if (args[n].dat == args[m].dat)
                already_done = true;
        }

        if (!already_done && args[n].opt && args[n].argtype == OP_ARG_DAT && args[n].dat->dirty_hd == Dirty::Host && 
            (!args[n].dat->set->is_particle || opp_current_device == Device_CPU)) 
        {
            oppic_download_dat(args[n].dat);
        }
    }

#ifdef USE_MPI
    return opp_mpi_halo_exchanges(set, nargs, args); // Halo exchange only for mesh dats
#else
    return set->size;
#endif
}

/*******************************************************************************
 * Routine to complete exchanging halos for the args
*******************************************************************************/
void opp_mpi_halo_wait_all(int nargs, oppic_arg *args)
{
#ifdef USE_MPI
    // Finish the halo exchange
    __opp_mpi_host_halo_wait_all(nargs, args); // Halo exchange only for mesh dats
#endif

    for (int n = 0; n < nargs; n++)
    { 
        bool already_done = false;

        // Check if dat reduction was already done within these args
        for (int m = 0; m < n; m++) 
        {
            if (args[n].dat == args[m].dat)
                already_done = true;
        }

        if (!already_done && args[n].opt && args[n].argtype == OP_ARG_DAT && opp_current_device == Device_GPU &&
                (args[n].dat->dirty_hd == Dirty::Device || !args[n].dat->set->is_particle)) 
        { 
            oppic_upload_dat(args[n].dat); // TODO : ideally need to upload only the halo regions               
        }
    }        

    for (int n = 0; n < nargs; n++) // This is wrong - change
    { 
        if (args[n].opt && args[n].argtype == OP_ARG_DAT && !args[n].dat->set->is_particle) 
        {
            if (opp_current_device == Device_CPU) 
                args[n].dat->dirty_hd = Dirty::Device;
            else
                args[n].dat->dirty_hd = Dirty::NotDirty;
        }
    }
}

