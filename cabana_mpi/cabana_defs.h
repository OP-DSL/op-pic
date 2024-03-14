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

#include <oppic_lib.h>
#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

#include "cabana_deck.h"

#define NG                        1
#define ONE                       1
#define TWO                       2
#define DIM                       3
#define ACCUMULATOR_ARRAY_LENGTH  4
#define INTERP_LEN                18
#define ACC_LEN (ACCUMULATOR_ARRAY_LENGTH*DIM)

#define VOXEL(x,y,z, nx,ny,nz) ((x) + ((nx)+(NG*2))*((y) + ((ny)+(NG*2))*(z)))

#define VOXEL_MAP(_ix,_iy,_iz, nx,ny,nz, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    int _z = _iz; \
    if ((_x < 0) || (_y < 0) || (_z < 0) || (_x > nx+NG) || (_y > ny+NG) || (_z > nz+NG)) OUT=-1; \
    else OUT = ((_x) + ((nx)+(NG*2))*((_y) + ((ny)+(NG*2))*(_z))); } \

#define VOXEL_MAP_NON_GHOST_PERIODIC(_ix,_iy,_iz, nx,ny,nz, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    int _z = _iz; \
    if (_x == 0) _x = (nx); \
    if (_y == 0) _y = (ny); \
    if (_z == 0) _z = (nz); \
    if (_x == nx+1) _x = 1; \
    if (_y == ny+1) _y = 1; \
    if (_z == nz+1) _z = 1; \
    OUT = ((_x) + ((nx)+(NG*2))*((_y) + ((ny)+(NG*2))*(_z))); } \

#define RANK_TO_INDEX(rank,ix,iy,iz,_x,_y) \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(_x);   /* ix = ix */                   \
    _iz  = _iy/int(_y);   /* iz = iz */                   \
    _iy -= _iz*int(_y);   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \

enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

enum CellAcc {
    jfx = 0 * ACCUMULATOR_ARRAY_LENGTH,
    jfy = 1 * ACCUMULATOR_ARRAY_LENGTH,
    jfz = 2 * ACCUMULATOR_ARRAY_LENGTH,
};

enum CellInterp {
    ex = 0,
    dexdy,
    dexdz,
    d2exdydz,
    ey,
    deydz,
    deydx,
    d2eydzdx,
    ez,
    dezdx,
    dezdy,
    d2ezdxdy,
    cbx,
    dcbxdx,
    cby,
    dcbydy,
    cbz,
    dcbzdz,
};

#define NEIGHBOURS 12
#define FACES 6

enum CellMap {
    xd_yd_z = 0,    // USED
    xd_y_zd,        // USED
    xd_y_z,         // USED
    x_yd_zd,        // USED
    x_yd_z,         // USED
    x_y_zd,         // USED
    x_y_zu,         // USED
    x_yu_z,         // USED
    x_yu_zu,        // USED
    xu_y_z,         // USED
    xu_y_zu,        // USED
    xu_yu_z,        // USED
};

enum Face {
    xd = 0,   
    yd = 1,  
    zd = 2,
    xu = 3,
    yu = 4,
    zu = 5,
};
