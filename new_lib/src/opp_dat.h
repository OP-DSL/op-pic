#pragma once

#include "opp_set.h"
#include "opp_arg.h"
#include "opp_util.h"

namespace opp
{
    class BaseDat
    {
    public:
        BaseDat() { };
        virtual ~BaseDat() { };
    };

    //*************************************************************************************************

    template<typename T> 
    class Dat : public BaseDat, public std::enable_shared_from_this<Dat<T>>
    {
    public:
        friend class Set;

        Dat(std::shared_ptr<Set> set, std::string name, UINT dimension, T* data);
        
        virtual ~Dat();
        virtual void print();

        inline UINT getIndex() { return dat_index; }
        inline std::shared_ptr<Set> getSet() { return set; }
        inline UINT getDimension() { return dimension; }
        inline UINT getSize() { return size; }
        inline std::vector<T>& getData() { return data; }
        inline std::string getDataType() { return type; }
        inline std::string getName() { return name; }
        inline opp_dirty getDirtyBit() { return dirtybit; }
        inline opp_dirty getDirtyHD() { return dirty_hd; }
        inline bool isUserManaged() { return user_managed; }

        inline void setIndex(UINT index) { dat_index = index; }
        
        inline T& get(INT m, INT n = 0) { return data[m * dimension + n]; }
        inline typename std::vector<T>::iterator& getIter() { return data.begin(); }

        Arg createArg(INT idx, std::shared_ptr<Map> map, INT dim, std::string type, opp_accessor acc, bool map_with_cell_index = false);
        Arg createArg(INT idx, std::shared_ptr<Map> map, opp_accessor acc, bool map_with_cell_index = false);
        Arg createArg(opp_accessor acc, bool map_with_cell_index = false);

    protected:
        UINT dat_index;
        std::shared_ptr<Set> set;
        UINT dimension;             /* dimension of data */
        UINT size;                  /* size of each element in dataset */
        std::vector<T> data;        /* data on host */
        T *data_d;                  /* data on device (GPU) */
        std::string type;           /* datatype */
        std::string name;           /* name of dataset */
        T *buffer_d;                /* buffer for MPI halo sends on the devidce */
        T *buffer_d_r;              /* buffer for MPI halo receives on the devidce */
        opp_dirty dirtybit;         /* flag to indicate MPI halo exchange is needed*/
        opp_dirty dirty_hd;         /* flag to indicate dirty status on host and device */
        bool user_managed;          /* indicates whether the user is managing memory */
        T *mpi_buffer;              /* ponter to hold the mpi buffer struct for the op_dat*/
    };

//*************************************************************************************************

    template<typename T> 
    class MeshDat : public virtual Dat<T> 
    {
        friend class MeshSet;

    public:
        virtual ~MeshDat();
        virtual void print();
    
    private:
        MeshDat(std::shared_ptr<MeshSet> set, std::string name, UINT dimension, T* data);

    };

    //*************************************************************************************************

    template<typename T> 
    class ParticleDat : public virtual Dat<T> 
    {
        friend class ParticleSet;

    public:  
        virtual ~ParticleDat();
        virtual void print();

    private:
        ParticleDat(std::shared_ptr<ParticleSet> set, std::string name, UINT dimension, opp_dat_type type, T* data);

        std::vector<std::vector<T>> thread_data;
        opp_dat_type dat_type;

    };

    //*************************************************************************************************

};