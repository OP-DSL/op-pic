#pragma once

#include "opp_types.h"
#include "opp_arg.h"

namespace opp
{
    class Map : public std::enable_shared_from_this<Map>
    {
        friend class MeshSet;

    public:
        virtual ~Map();

        void print();

        inline UINT getIndex() { return map_index; }
        inline UINT getDimension() { return dimension; }
        inline std::shared_ptr<Set> getFrom() { return from_set; }
        inline std::shared_ptr<Set> getTo() { return to_set; }
        inline const std::vector<INT>& getMapData() { return map; }
        inline std::string getName() { return name; }
        inline INT get(INT m, INT n = 0) { return map[m * dimension + n]; }

        Arg createArg(opp_accessor acc, bool map_with_cell_index = false);

    protected:
        Map(std::string name, std::shared_ptr<Set> from_set, std::shared_ptr<Set> to_set, INT dimension, INT* map);
        
        UINT map_index;
        UINT dimension;
        std::string name;
        std::shared_ptr<Set> from_set;
        std::shared_ptr<Set> to_set;
        std::vector<INT> map;
        INT *map_d = nullptr;
    };
};