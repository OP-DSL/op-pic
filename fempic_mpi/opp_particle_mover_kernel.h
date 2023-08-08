#pragma once

#include <oppic_lib.h>

const double ONE_OVER_SIX = (1.0 / 6.0);
const int N_PER_C = 4;
const int DET_FIELDS = 4;
const int NEIGHB_C = 4;

//*******************************************************************************
inline bool isPointInCellKernel(const double *point_pos, int* current_cell_index,
    double* point_lc, const double *cell_volume, const double *cell_det) { 

    bool inside = true;  
    double coefficient2 = ONE_OVER_SIX / (*cell_volume);

    for (int i=0; i<N_PER_C; i++) { /*loop over vertices*/
    
        point_lc[i] = coefficient2 * (
            cell_det[i * DET_FIELDS + 0] - 
            cell_det[i * DET_FIELDS + 1] * point_pos[0] + 
            cell_det[i * DET_FIELDS + 2] * point_pos[1] - 
            cell_det[i * DET_FIELDS + 3] * point_pos[2]);
        
        if (point_lc[i] < 0.0 || 
            point_lc[i] > 1.0)  
                inside = false;
    }    
    
    return inside;
}

//*******************************************************************************
inline opp_move_status getCellIndexKernel(const double *point_pos, int* current_cell_index,
    double* point_lc, const double *cell_volume, const double *cell_det, const int *cell_connectivity) { 

    if (isPointInCellKernel(point_pos, current_cell_index, point_lc, cell_volume, cell_det)) {
        return OPP_MOVE_DONE;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = point_lc[0];
    
    for (int i=1; i<NEIGHB_C; i++) {
        if (point_lc[i] < min_lc) {
            min_lc = point_lc[i];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i] >= 0) { // is there a neighbor in this direction?
        (*current_cell_index) = cell_connectivity[min_i];
        return OPP_NEED_MOVE;
    }
    else {
        (*current_cell_index) = MAX_CELL_INDEX;
        return OPP_NEED_REMOVE;
    }
}