#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

#define INT int
#define UINT unsigned int
#define REAL double
#define BOOL bool

#define ZERO_double    0.0
#define ZERO_float     0.0f
#define ZERO_int       0
#define ZERO_uint      0
#define ZERO_ll        0
#define ZERO_ull       0
#define ZERO_bool      0

#define INT_IT std::vector<INT>::iterator
#define REAL_IT std::vector<REAL>::iterator

#define Ptr(a) std::shared_ptr<a>
#define Cast_X_to_Type(X, Type) std::dynamic_pointer_cast<Type>(X)

namespace opp
{
    enum opp_accessor
    {
        OPP_READ = 0,
        OPP_WRITE,
        OPP_RW,   
        OPP_INC, 
        OPP_MIN,  
        OPP_MAX,  
    };

    enum opp_arg_type
    {
        OPP_ARG_GBL = 0,
        OPP_ARG_DAT,
        OPP_ARG_MAP,   
    };

    enum opp_dat_type
    {
        OPP_DAT_GENERAL = 0,
        OPP_DAT_CELL_INDEX,
        OPP_DAT_POSITION,   
    };

    enum opp_iterate_type
    {
        OPP_ITERATE_ALL = 1,
        OPP_ITERATE_INJ,
    };

    enum opp_move_status 
    {
        OPP_MOVE_DONE = 0,
        OPP_NEED_MOVE,
        OPP_NEED_REMOVE,
    };

    enum opp_device_type
    {
        OPP_DEVICE_CPU = 1,
        OPP_DEVICE_GPU = 2,
    };

    enum opp_dirty
    {
        OPP_NOT_DIRTY = 0,
        OPP_DEVICE_DIRTY = 1,
        OPP_HOST_DIRTY = 2,
    };

    class BaseDat;
    template<typename T> class Dat;
    template<typename T> class ParticleDat;
    template<typename T> class MeshDat;
    class Map;
    class Set;
    class MeshSet;
    class ParticleSet;
    class Simulation;

    struct opp_arg 
    {
        std::shared_ptr<BaseDat> dat;   /* dataset */
        std::shared_ptr<Map> map;       /* indirect mapping */
        char const *type;               /* datatype */
        opp_accessor acc;               /* opp_accessor OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
        opp_arg_type argtype;           /* argument type - global, map, dat */
        int sent;                       /* flag to indicate if this argument has data in flight under non-blocking MPI comms*/
    };
};