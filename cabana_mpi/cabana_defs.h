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
    #include <opp_mpi.h>
#endif

#ifdef DEBUG_LOG
    #define FP_DEBUG true
#else
    #define FP_DEBUG false
#endif

#define ONE                       1
#define TWO                       2
#define DIM                       3
#define ACCUMULATOR_ARRAY_LENGTH  4
#define INTERP_LEN                18

#define VOXEL(x,y,z, nx,ny,nz, NG) ((x) + ((nx)+(NG*2))*((y) + ((ny)+(NG*2))*(z)))

#define VOXEL_MAP(_ix,_iy,_iz, nx,ny,nz, NG, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    int _z = _iz; \
    if (_x < 0) _x = (nx+2*NG-1); \
    if (_y < 0) _y = (ny+2*NG-1); \
    if (_z < 0) _z = (nz+2*NG-1); \
    if (_x > nx+NG) _x = 0; \
    if (_y > ny+NG) _y = 0; \
    if (_z > nz+NG) _z = 0; \
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

#define NEIGHBOUR_CELLS 26

enum CellMap {
    xd_yd_zd = 0,
    xd_yd_z,
    xd_yd_zu,
    xd_y_zd,
    xd_y_z,
    xd_y_zu,
    xd_yu_zd,
    xd_yu_z,
    xd_yu_zu,
    x_yd_zd,
    x_yd_z,
    x_yd_zu,
    x_y_zd,
    x_y_zu,
    x_yu_zd,
    x_yu_z,
    x_yu_zu,
    xu_yd_zd,
    xu_yd_z,
    xu_yd_zu,
    xu_y_zd,
    xu_y_z,
    xu_y_zu,
    xu_yu_zd,
    xu_yu_z,
    xu_yu_zu,
    interp_0,
    interp_1,
    interp_2,
    interp_3,
    interp_4,
    interp_5,
};

enum FACE {
    FACE_X_MIN  = 0,
    FACE_X_PLUS = 1,
    FACE_Y_MIN  = 2,
    FACE_Y_PLUS = 3,
    FACE_Z_MIN  = 4,
    FACE_Z_PLUS = 5,
};


//*************************************************************************************************
inline int detect_leaving_domain(size_t ix, size_t iy, size_t iz, int nx, int ny, int nz) {

    int leaving = -1;
    if (ix == 0) leaving = 0;
    if (iy == 0) leaving = 1;
    if (iz == 0) leaving = 2;
    if (ix == nx+1) leaving = 3;
    if (iy == ny+1) leaving = 4;
    if (iz == nz+1) leaving = 5;

    return leaving;
}

//*************************************************************************************************
inline int get_neighbour_cell(int cell_index, FACE face, int nx, int ny, int nz, int ng)
{
	size_t ix, iy, iz;
	RANK_TO_INDEX(cell_index, ix, iy, iz, (nx+(2*ng)), (ny+(2*ng)));

	if (face == FACE_X_MIN ) { ix--; }
	else if (face == FACE_X_PLUS) { iy--; }
	else if (face == FACE_Y_MIN ) { iz--; }
	else if (face == FACE_Y_PLUS) { ix++; }
	else if (face == FACE_Z_MIN ) { iy++; }
	else if (face == FACE_Z_PLUS) { iz++; }

	const int is_leaving_domain = detect_leaving_domain(ix, iy, iz, nx, ny, nz);
	if (is_leaving_domain >= 0) {

		//if ( boundary == Boundary::Periodic)
		{
			if (is_leaving_domain == 0) { // -1 on x face
				ix = (nx-1) + ng;
			}
			else if (is_leaving_domain == 1) { // -1 on y face
				iy = (ny-1) + ng;
			}
			else if (is_leaving_domain == 2) { // -1 on z face
				iz = (nz-1) + ng;
			}
			else if (is_leaving_domain == 3) { // 1 on x face
				ix = ng;
			}
			else if (is_leaving_domain == 4) { // 1 on y face
				iy = ng;
			}
			else if (is_leaving_domain == 5) { // 1 on z face
				iz = ng;
			}
		}
	}

	return VOXEL(ix, iy, iz, nx, ny, nz, ng);
}


void opp_loop_all__interpolate_mesh_fields(
    opp_set set,        // cells_set
    opp_arg arg0,       // cell0_e,        // OPP_READ
    opp_arg arg1,       // cell0_b,        // OPP_READ
    opp_arg arg2,       // cell0_ghost,    // OPP_READ
    opp_arg arg3,       // cell_x_e,       // OPP_READ
    opp_arg arg4,       // cell_y_e,       // OPP_READ
    opp_arg arg5,       // cell_z_e,       // OPP_READ
    opp_arg arg6,       // cell_yz_e,      // OPP_READ 
    opp_arg arg7,       // cell_xz_e,      // OPP_READ
    opp_arg arg8,       // cell_xy_e,      // OPP_READ
    opp_arg arg9,       // cell_x_b,       // OPP_READ
    opp_arg arg10,      // cell_y_b,       // OPP_READ
    opp_arg arg11,      // cell_z_b        // OPP_READ
    opp_arg arg12       // cell0_interp    // OPP_WRITE
);

void opp_particle_mover__Move(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_cid           // OPP_RW
    opp_arg arg1,       // part_vel           // OPP_RW
    opp_arg arg2,       // part_pos           // OPP_RW
    opp_arg arg3,       // part_streak_mid    // OPP_RW
    opp_arg arg4,       // part_weight        // OPP_READ
    opp_arg arg5,       // cell_inter         // OPP_READ
    opp_arg arg6        // cell_acc           // OPP_INC
);

void opp_loop_all__accumulate_current_to_cells(
    opp_set set,     // cells set
    opp_arg arg0,    // iter_acc        // OPP_READ
    opp_arg arg1,    // cell0_j         // OPP_WRITE
    opp_arg arg2,    // cell0_acc       // OPP_READ
    opp_arg arg3,    // cell_xd_acc     // OPP_READ
    opp_arg arg4,    // cell_yd_acc     // OPP_READ
    opp_arg arg5,    // cell_zd_acc     // OPP_READ
    opp_arg arg6,    // cell_xyd_acc    // OPP_READ
    opp_arg arg7,    // cell_yzd_acc    // OPP_READ
    opp_arg arg8     // cell_xzd_acc    // OPP_READ
);

void opp_loop_all__half_advance_b(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_ghost     // OPP_READ
    opp_arg arg1,    // cell_x_e        // OPP_WRITE
    opp_arg arg2,    // cell_y_e        // OPP_READ
    opp_arg arg3,    // cell_z_e        // OPP_READ
    opp_arg arg4,    // cell0_e         // OPP_READ
    opp_arg arg5     // cell0_b         // OPP_INC
);

void opp_loop_all__advance_e(
    opp_set set,     // cells set
    opp_arg arg0,    // iter_adv_e      // OPP_READ
    opp_arg arg1,    // cell_x_b        // OPP_READ
    opp_arg arg2,    // cell_y_b        // OPP_READ
    opp_arg arg3,    // cell_z_b        // OPP_READ
    opp_arg arg4,    // cell0_b         // OPP_READ
    opp_arg arg5,    // cell0_j         // OPP_READ
    opp_arg arg6     // cell0_e         // OPP_INC
);
