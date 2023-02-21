#pragma once

#include "opp_types.h"

namespace opp
{
    class Arg 
    {
    public:
        Arg(std::shared_ptr<BaseDat> dat, INT idx, std::shared_ptr<Map> map, INT dim, std::string type, opp_accessor acc, opp_arg_type argtype, bool map_with_cell_index = false);
        virtual ~Arg();

        inline const std::vector<INT>& getMapData() { return *(i_map); }
        inline std::vector<INT>& getIData() { return *i_dat; }
        inline std::vector<REAL>& getRData() { return *r_dat; }
        inline std::vector<BOOL>& getBData() { return *b_dat; }

        // DO I NEED TO REMOVE THE BELOW?
        inline INT& datI(INT m, INT n = 0) { return (*i_dat)[m * dim + n]; }
        inline REAL& datR(INT m, INT n = 0) { return (*r_dat)[m * dim + n]; }
        // inline BOOL& datB(INT m, INT n = 0) { return (*b_dat)[m * dim + n]; }

        std::shared_ptr<BaseDat> dat;                   /* dataset */
        std::shared_ptr<Map> map;                       /* indirect mapping */
        std::string type            = "";               /* datatype */
        INT dim                     = 0;
        INT idx                     = 0;
        bool map_with_cell_index    = false;
        opp_accessor acc            = OPP_READ;         /* opp_accessor OP_READ, OP_WRITE, OP_RW, OP_INC, OP_MIN, OP_MAX */
        opp_arg_type argtype        = OPP_ARG_DAT;      /* argument type - global, map, dat */
        int sent                    = 0;                /* flag to indicate if this argument has data in flight under non-blocking MPI comms*/

    protected:
        const std::vector<INT>* i_map = nullptr;
        std::vector<INT>* i_dat       = nullptr;
        std::vector<REAL>* r_dat      = nullptr;
        std::vector<BOOL>* b_dat      = nullptr;
    };
};
