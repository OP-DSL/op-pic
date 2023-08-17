#pragma once

#include <oppic_lib.h>

const double ONE_OVER_SIX = (1.0 / 6.0);
const int N_PER_C = 4;
const int DET_FIELDS = 4;
const int NEIGHB_C = 4;

#define OPP_LOOP_UNROLL
//*******************************************************************************
inline bool isPointInCellKernel(const double *point_pos, double* point_lc, 
                                const double *cell_volume, const double *cell_det) { 

    bool inside = true;  
    double coefficient2 = ONE_OVER_SIX / (*cell_volume);

#ifdef OPP_LOOP_UNROLL

    // alignas(32) double pos0 = point_pos[0];
    // alignas(32) double pos1 = point_pos[1];
    // alignas(32) double pos2 = point_pos[2];

    // alignas(32) double det0 = cell_det[0];
    // alignas(32) double det1 = cell_det[1];
    // alignas(32) double det2 = cell_det[2];
    // alignas(32) double det3 = cell_det[3];
    // alignas(32) double det4 = cell_det[4];
    // alignas(32) double det5 = cell_det[5];
    // alignas(32) double det6 = cell_det[6];
    // alignas(32) double det7 = cell_det[7];
    // alignas(32) double det8 = cell_det[8];
    // alignas(32) double det9 = cell_det[9];
    // alignas(32) double det10 = cell_det[10];
    // alignas(32) double det11 = cell_det[11];
    // alignas(32) double det12 = cell_det[12];
    // alignas(32) double det13 = cell_det[13];
    // alignas(32) double det14 = cell_det[14];
    // alignas(32) double det15 = cell_det[15];

    // alignas(32) double point_lc0 = coefficient2 * (
    //         det0 - 
    //         det1 * pos0 + 
    //         det2 * pos1 - 
    //         det3 * pos2);

    // alignas(32) double point_lc1 = coefficient2 * (
    //         det4 - 
    //         det5 * pos0 + 
    //         det6 * pos1 - 
    //         det7 * pos2);

    // alignas(32) double point_lc2 = coefficient2 * (
    //         det8  - 
    //         det9  * pos0 + 
    //         det10 * pos1 - 
    //         det11 * pos2);

    // alignas(32) double point_lc3 = coefficient2 * (
    //         det12 - 
    //         det13 * pos0 + 
    //         det14 * pos1 - 
    //         det15 * pos2);

    // if (point_lc0 < 0.0 || 
    //     point_lc0 > 1.0 ||
    //     point_lc1 < 0.0 || 
    //     point_lc1 > 1.0 ||
    //     point_lc2 < 0.0 || 
    //     point_lc2 > 1.0 ||
    //     point_lc3 < 0.0 || 
    //     point_lc3 > 1.0) { 
            
    //         inside = false;
    // }

    // point_lc[0] = point_lc0;
    // point_lc[1] = point_lc1;
    // point_lc[2] = point_lc2;
    // point_lc[3] = point_lc3;

    point_lc[0] = coefficient2 * (
            cell_det[0] - 
            cell_det[1] * point_pos[0] + 
            cell_det[2] * point_pos[1] - 
            cell_det[3] * point_pos[2]);

    point_lc[1] = coefficient2 * (
            cell_det[4] - 
            cell_det[5] * point_pos[0] + 
            cell_det[6] * point_pos[1] - 
            cell_det[7] * point_pos[2]);

    point_lc[2] = coefficient2 * (
            cell_det[8] - 
            cell_det[9] * point_pos[0] + 
            cell_det[10] * point_pos[1] - 
            cell_det[11] * point_pos[2]);

    point_lc[3] = coefficient2 * (
            cell_det[12] - 
            cell_det[13] * point_pos[0] + 
            cell_det[14] * point_pos[1] - 
            cell_det[15] * point_pos[2]);

    if (point_lc[0] < 0.0 || 
        point_lc[0] > 1.0 ||
        point_lc[1] < 0.0 || 
        point_lc[1] > 1.0 ||
        point_lc[2] < 0.0 || 
        point_lc[2] > 1.0 ||
        point_lc[3] < 0.0 || 
        point_lc[3] > 1.0) { 
            
            inside = false;
    }

#else
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
#endif
    
    return inside;
}

//*******************************************************************************
inline opp_move_status getCellIndexKernel(const double *point_pos, int* current_cell_index,
    double* point_lc, const double *cell_volume, const double *cell_det, const int *cell_connectivity) { 

    if (isPointInCellKernel(point_pos, point_lc, cell_volume, cell_det)) {
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