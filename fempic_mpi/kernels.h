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

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************


#include "fempic_ori/maths.h"
#include "oppic_lib.h"

const double K_ONE_OVER_SIX = (1.0 / 6.0);

//*************************************************************************************************
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

//*************************************************************************************************
inline void inject_ions__kernel(
    double *part_pos,
    double *part_vel,
    int *part_cell_connectivity,
    int *cell_id, 
    double *cell_ef,
    double *iface_u,
    double *iface_v,
    double *iface_normal,
    double *node_pos,
    const double* dummy_part_random
)
{
    double a = dummy_part_random[0];
    double b = dummy_part_random[1];
    if ((a + b) > 1)  
    {
        a = (1 - a);
        b = (1 - b);
    }

    for (int i = 0; i < DIM; i++) 
    {
        part_pos[i] = a * iface_u[i] + b * iface_v[i] + node_pos[i];
        
        part_vel[i] = (iface_normal[i] * CONST_ion_velocity);
        part_vel[i] -= CONST_charge / CONST_mass * cell_ef[i] * (0.5 * CONST_dt);
    }

    (*part_cell_connectivity) = (*cell_id);
}

inline void calculate_new_pos_vel__kernel(
    const double *cell_ef,
    double *part_pos,
    double *part_vel ) {

    double coefficient1 = CONST_charge / CONST_mass * (CONST_dt);
    for (int i = 0; i < DIM; i++)
        part_vel[i] += (coefficient1 * cell_ef[i]);                  
    
    for (int i = 0; i < DIM; i++)
        part_pos[i] += part_vel[i] * (CONST_dt); // v = u + at
}

inline void deposit_charge_on_nodes__kernel(
    const double *part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3) {

    (*node_charge_den0) += part_lc[0];
    (*node_charge_den1) += part_lc[1];
    (*node_charge_den2) += part_lc[2];
    (*node_charge_den3) += part_lc[3];
}

//*************************************************************************************************
inline void move_all_particles_to_cell__kernel(
    opp_move_var& m,
    const double *cell_ef,
    double *part_pos,
    double *part_vel,
    double *part_lc,
    int* current_cell_index,
    const double *current_cell_volume,
    const double *current_cell_det,
    const int *cell_connectivity,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
)
{
    if (m.iteration_one)
    {
        double coefficient1 = CONST_charge / CONST_mass * (CONST_dt);
        for (int i = 0; i < DIM; i++)
            part_vel[i] += (coefficient1 * cell_ef[i]);                  
        
        for (int i = 0; i < DIM; i++)
            part_pos[i] += part_vel[i] * (CONST_dt); // v = u + at
    }

    bool inside = true;
    double coefficient2 = K_ONE_OVER_SIX / (*current_cell_volume);
    for (int i=0; i<N_PER_C; i++) /*loop over vertices*/
    {
        part_lc[i] = coefficient2 * (
            current_cell_det[i * DET_FIELDS + 0] - 
            current_cell_det[i * DET_FIELDS + 1] * part_pos[0] + 
            current_cell_det[i * DET_FIELDS + 2] * part_pos[1] - 
            current_cell_det[i * DET_FIELDS + 3] * part_pos[2]);
        
        if (part_lc[i] < 0.0 || 
            part_lc[i] > 1.0)  
                inside = false;
                // m.inside_cell = false;
    }    
    
    // if (m.inside_cell)
    if (inside)
    {
        m.move_status = OPP_MOVE_DONE;

        (*node_charge_den0) += part_lc[0];
        (*node_charge_den1) += part_lc[1];
        (*node_charge_den2) += part_lc[2];
        (*node_charge_den3) += part_lc[3];

        return;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0];
    
    for (int i=1; i<NEIGHB_C; i++)
    {
        if (part_lc[i] < min_lc) 
        {
            min_lc = part_lc[i];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i] >= 0) // is there a neighbor in this direction?
    {
        (*current_cell_index) = cell_connectivity[min_i];
        m.move_status = OPP_NEED_MOVE;
    }
    else
    {
        m.move_status = OPP_NEED_REMOVE;
    }
}

//*************************************************************************************************
inline void compute_node_charge_density__kernel(
    double *node_charge_den,
    const double *node_volume
)
{
    (*node_charge_den) *= (CONST_spwt / (*node_volume));
}

//*************************************************************************************************
inline void compute_electric_field__kernel(
    double *cell_electric_field,             
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
)
{
    for (int dim = 0; dim < DIM; dim++)
    { 
        double c1 = (cell_shape_deriv[0 * DIM + dim] * (*node_potential0));
        double c2 = (cell_shape_deriv[1 * DIM + dim] * (*node_potential1));
        double c3 = (cell_shape_deriv[2 * DIM + dim] * (*node_potential2));
        double c4 = (cell_shape_deriv[3 * DIM + dim] * (*node_potential3));

        cell_electric_field[dim] -= (c1 + c2 + c3 + c4);
    }    
}

//*************************************************************************************************