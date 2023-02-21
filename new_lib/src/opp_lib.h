#pragma once

#include "opp_set.h"
#include "opp_map.h"
#include "opp_dat.h"
#include "opp_arg.h"

#include "opp_simulation.h"

namespace opp
{
    template<typename T> 
    std::shared_ptr<MeshDat<T>> MeshSet::createDat(std::string name, INT dimension, T* data)
    {
        std::shared_ptr<MeshDat<T>> meshDat(new MeshDat<T>(this->shared_from_this(), std::move(name), dimension, data));
        meshDat->setIndex(SIM->registerDat(meshDat)); 
        this->dats.push_back(meshDat);
        
        return meshDat;
    }

    template<typename T> 
    std::shared_ptr<ParticleDat<T>> ParticleSet::createDat(std::string name, INT dimension, T* data, opp_dat_type type)
    {
        std::shared_ptr<ParticleDat<T>> particleDat(new ParticleDat<T>(this->shared_from_this(), std::move(name), dimension, type, data));

        if (type == OPP_DAT_CELL_INDEX)
        {
            if constexpr (std::is_same<T, INT>::value)
            {
                this->setCellIndexDat(particleDat);
            }
            else
            {
                printf("ParticleSet::createDat error -- cell index provided %s should be of type INT\n", name.c_str());
                exit(-1);                    
            }
        }

        return particleDat;
    }

    template<typename T> 
    Dat<T>::Dat(std::shared_ptr<Set> set, std::string name, UINT dimension, T* data)
            : set(set), name(std::move(name)), dimension(dimension)
    {

        if (dimension == 0) {
            printf("Set::enrichDat error -- zero dimension for data: %s\n", name.c_str());
            exit(-1);
        }

        if constexpr (std::is_same<T, INT>::value)
        {
            this->type = "int";
            this->size = dimension * sizeof(INT);
        }
        else if constexpr (std::is_same<T, REAL>::value)
        {
            this->type = "double";
            this->size = dimension * sizeof(REAL);
        }
        else
        {
            printf("Set::enrichDat error -- support for data type not implemented for data: %s\n", name.c_str());
            exit(-1);
        }

        this->user_managed  = false;
        this->mpi_buffer    = nullptr;
        this->buffer_d      = nullptr;
        this->buffer_d_r    = nullptr;
        this->dirty_hd      = OPP_NOT_DIRTY;
        this->dirtybit      = OPP_NOT_DIRTY;

        if (data != nullptr)
            this->data = util::getVecFromPtr<T>(data, dimension * set->getSize());
        else if (set->getSize() > 0)
            this->data.reserve(dimension * set->getSize());
    }

    template<typename T> Dat<T>::~Dat()
    {
        std::cout << "At distructor of [Dat] " << name << std::endl;
    }

    template<typename T> void Dat<T>::print() 
    {
        std::cout << name << " [Dat] size: " << size << std::endl; 
    }

    template<typename T> Arg Dat<T>::createArg(INT idx, std::shared_ptr<Map> map, INT dim, std::string type, opp_accessor acc, bool map_with_cell_index)
    {
        return Arg(this->shared_from_this(), idx, map, dim, type, acc, OPP_ARG_DAT, map_with_cell_index);
    }

    template<typename T> Arg Dat<T>::createArg(INT idx, std::shared_ptr<Map> map, opp_accessor acc, bool map_with_cell_index)
    {
        return Arg(this->shared_from_this(), idx, map, this->dimension, this->type, acc, OPP_ARG_DAT, map_with_cell_index);
    }

    template<typename T> Arg Dat<T>::createArg(opp_accessor acc, bool map_with_cell_index)
    {
        return Arg(this->shared_from_this(), 0, std::shared_ptr<Map>(nullptr), this->dimension, this->type, acc, OPP_ARG_DAT, map_with_cell_index);
    }

    template<typename T> 
    ParticleDat<T>::ParticleDat(std::shared_ptr<ParticleSet> set, std::string name, UINT dimension, opp_dat_type type, T* data)
            : Dat<T>(std::static_pointer_cast<Set>(set), name, dimension, data)
    {
        std::cout << "At constructor of [ParticleDat] " << this->name << " size: " << this->size << std::endl;
        
        this->dat_type = type;
    }

    template<typename T> ParticleDat<T>::~ParticleDat()
    {
        std::cout << "At distructor of [ParticleDat] " << this->name << std::endl;
    }

    template<typename T> void ParticleDat<T>::print() 
    {
        std::cout << this->name << " [ParticleDat] size: " << this->size << std::endl; 
    }

    template<typename T> MeshDat<T>::MeshDat(std::shared_ptr<MeshSet> set, std::string name, UINT dimension, T* data)
            : Dat<T>(std::static_pointer_cast<Set>(set), name, dimension, data)
    {
        std::cout << "At constructor of [MeshDat] " << this->name << " size: " << this->size << std::endl;
    }

    template<typename T> MeshDat<T>::~MeshDat()
    {
        std::cout << "At distructor of [MeshDat] " << this->name << std::endl;
    }

    template<typename T> void MeshDat<T>::print() 
    {
        std::cout << this->name << " [MeshDat] size: " << this->size << std::endl; 
    }
};