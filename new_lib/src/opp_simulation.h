#pragma once

#include "opp_set.h"
#include "opp_map.h"
#include "opp_dat.h"
#include "opp_arg.h"

#define SIM Simulation::getInstance()

namespace opp
{
    class Simulation
    {
    public:
        Simulation(const Simulation& obj) = delete;

        static Simulation* getInstance();
        static void releaseInstance();

        void setValues(std::string name, std::string loves);
        void print();
    
        int OPP_maps_base_index = 0;

        // template<typename U> 
        // inline Arg createGlobalArg(U* data, INT dim, opp_accessor acc)
        // {
        //     return Arg(this->shared_from_this(), idx, map, dim, type, acc, OPP_ARG_DAT, map_with_cell_index);
        // }
        
        std::shared_ptr<MeshSet> createMeshSet(std::string name, int size);
        std::shared_ptr<ParticleSet> createParticleSet(std::string name, std::shared_ptr<MeshSet> cells_set, int size = 0);

        UINT registerDat(std::shared_ptr<BaseDat> dat);
        UINT registerMap(std::shared_ptr<Map> map);
        UINT registerSet(std::shared_ptr<Set> set);

    private:

        Simulation();   
        virtual ~Simulation();

        static Simulation* instancePtr;

        // TODO : enrich below and index of set, map, dat
        std::vector<std::shared_ptr<Set>> opp_sets;
        std::vector<std::shared_ptr<Map>> opp_maps;
        std::vector<std::shared_ptr<BaseDat>> opp_dats;

        std::string name, loves;
    };

};