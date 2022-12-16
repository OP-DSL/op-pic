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

#include "fempic.h"

//*************************************************************************************************
void inject_ions__kernel(
    double *pos,
    double *vel,
    double *ef ,
    double *lc, 
    int *cell_index
)
{
    /*sample random position on the inlet face*/
    pos[0] = -0.1 + 0.2 * rnd();
    pos[1] = -0.1 + 0.2 * rnd();
    pos[2] = 0;

    /*injecting cold beam*/
    vel[0] = 0;
    vel[1] = 0;
    vel[2] = ION_VELOCITY;

    ef[0] = 0;
    ef[1] = 0;
    ef[2] = 0;

    lc[0] = 0.0;
    lc[1] = 0.0;
    lc[2] = 0.0;
    lc[3] = 0.0;

    *cell_index = 0;
}

//*************************************************************************************************
void enrich_velocity__kernel(
    double *vel,
    const double *cell_ef,
    double *dt
)
{
    vel[0] -= OP_CONST_charge / OP_CONST_mass * cell_ef[0] * (0.5 * (*dt));
    vel[1] -= OP_CONST_charge / OP_CONST_mass * cell_ef[1] * (0.5 * (*dt));
    vel[2] -= OP_CONST_charge / OP_CONST_mass * cell_ef[2] * (0.5 * (*dt));
}

//*************************************************************************************************
void move_particle_to_cell__kernel(
    int* move_status,
    const double* part_pos,
    double* part_lc,
    int* current_cell_index,
    const double *cell_volume,
    const double *cell_det,
    const int *cell_connectivity,
    const bool* search)
{
    bool inside = true;

    for (int i=0; i<NODES_PER_CELL; i++) /*loop over vertices*/
    {
        part_lc[i] = (1.0/6.0) * (
            cell_det[i * DET_FIELDS + 0] - 
            cell_det[i * DET_FIELDS + 1] * part_pos[0] + 
            cell_det[i * DET_FIELDS + 2] * part_pos[1] - 
            cell_det[i * DET_FIELDS + 3] * part_pos[2]
                ) / (*cell_volume);
        
        if (part_lc[i]<0 || part_lc[i]>1.0) inside = false;
    }    
    
    if (inside)
    {
        *move_status = MOVE_DONE;
        return;
    }

    if (*search) 
    {
        (*current_cell_index)++; // outside the last known cell, Increment the cell_index to search in the full mesh
        return;
    }

    // outside the last known cell, find most negative weight and use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0];
    
    for (int i=1; i<NEIGHBOUR_CELLS; i++)
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
        *move_status = NEED_MOVE;
    }
    else
    {
        (*current_cell_index) = -1;
        *move_status = NEED_REMOVE;
    }
}

//*************************************************************************************************
void weight_fields_to_particles__kernel(
    double *part_ef,
    const double *cell_ef
)
{
    part_ef[0] = cell_ef[0];
    part_ef[1] = cell_ef[1];
    part_ef[2] = cell_ef[2];
}

//*************************************************************************************************
void move_particles__kernel(
    double *pos,    
    double *vel,
    const double *ef,
    const double *dt
)
{
    vel[0] += (OP_CONST_charge / OP_CONST_mass * ef[0] * (*dt));
    vel[1] += (OP_CONST_charge / OP_CONST_mass * ef[1] * (*dt));
    vel[2] += (OP_CONST_charge / OP_CONST_mass * ef[2] * (*dt));
    
    pos[0] += vel[0] * (*dt); // v = u + at
    pos[1] += vel[1] * (*dt); // v = u + at
    pos[2] += vel[2] * (*dt); // v = u + at
}

//*************************************************************************************************
void reset_ion_density__kernel(
    double *ion_den
)
{
    ion_den[0] = 0.0;
}

//*************************************************************************************************
void weight_particle_to_mesh_nodes__kernel(
    const double* part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3,
    const double *node_volume0,
    const double *node_volume1,
    const double *node_volume2,
    const double *node_volume3
)
{
    (*node_charge_den0) += (part_lc[0] * (OP_CONST_spwt / (*node_volume0)));
    (*node_charge_den1) += (part_lc[1] * (OP_CONST_spwt / (*node_volume1)));
    (*node_charge_den2) += (part_lc[2] * (OP_CONST_spwt / (*node_volume2)));
    (*node_charge_den3) += (part_lc[3] * (OP_CONST_spwt / (*node_volume3)));
}

//*************************************************************************************************