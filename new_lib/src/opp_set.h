#pragma once

#include "opp_types.h"
#include "opp_util.h"
#include "opp_map.h"

namespace opp
{
    class Set
    {
    public:
        virtual void print();

        virtual ~Set();

        inline UINT getIndex() { return set_index; }
        inline UINT getSize() { return size; }
        inline UINT getCapacity() { return capacity; }
        inline UINT getCoreSize() { return core_size; }
        inline UINT getExecSize() { return exec_size; }
        inline UINT getNonExecSize() { return nonexec_size; }
        inline std::string getName() { return name; }

    protected:
        Set(std::string name, int size = 0);

        std::string name  = "";
        UINT set_index    = 0;
        UINT size         = 0;
        UINT capacity     = 0;
        UINT core_size    = 0;
        UINT exec_size    = 0;
        UINT nonexec_size = 0;

        std::vector<std::shared_ptr<BaseDat>> dats;
        std::vector<std::shared_ptr<Map>> maps; 
    };

    //*************************************************************************************************

    class MeshSet : public Set, public std::enable_shared_from_this<MeshSet> 
    {
    public:
        friend class MeshDat<INT>;
        friend class MeshDat<REAL>;
        friend class Map;

        MeshSet(std::string name, int size);
        virtual ~MeshSet();

        virtual void print() override;

        template<typename T> 
        std::shared_ptr<MeshDat<T>> createDat(std::string name, INT dimension, T* data = nullptr);

        std::shared_ptr<Map> createMapTo(std::shared_ptr<MeshSet> to_set, std::string name, INT dim, INT* data);
    };

    //*************************************************************************************************

    class ParticleSet : public Set, public std::enable_shared_from_this<ParticleSet>
    {
    public:
        ParticleSet(std::string name, std::shared_ptr<MeshSet> cells_set, int size = 0);
        virtual ~ParticleSet();

        virtual void print() override;
        virtual void increaseParticleCount(UINT add_count);

        inline void resetAddCount() { add_count = 0; }
        inline std::shared_ptr<ParticleDat<INT>> getCellIndexDat() { return cell_index_dat; }
        inline std::shared_ptr<MeshSet> getCellsSet() { return cells_set; }
        inline void setCellIndexDat(std::shared_ptr<ParticleDat<INT>> dat) { cell_index_dat = dat; } 

        template<typename T> 
        std::shared_ptr<ParticleDat<T>> createDat(std::string name, INT dimension, T* data = nullptr, opp_dat_type type = OPP_DAT_GENERAL);

    protected:
        UINT add_count = 0;
        std::shared_ptr<ParticleDat<INT>> cell_index_dat;
        std::shared_ptr<MeshSet> cells_set;
    };
};