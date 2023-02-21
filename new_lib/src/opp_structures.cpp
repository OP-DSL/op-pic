
#include "opp_lib.h"

using namespace opp;

//*************************************************************************************************

Set::Set(std::string name, int size)
    : name(std::move(name)), size(size)
{
    std::cout << "At constructor of [Set] " << name << " size: " << size << std::endl;
}

Set::~Set()
{
    std::cout << "At distructor of [Set] " << name << std::endl;
}

void Set::print() 
{
    std::cout << name << " [Set] size: " << size << std::endl; 
}

//*************************************************************************************************

MeshSet::MeshSet(std::string name, int size)
    : Set(std::move(name), size)
{
    std::cout << "At constructor of [MeshSet] " << this->name << " size: " << this->size << std::endl;
}

MeshSet::~MeshSet()
{
    std::cout << "At distructor of [MeshSet] " << this->name << std::endl;
}

void MeshSet::print() 
{
    std::cout << this->name << " [MeshSet] size: " << this->size << std::endl; 
}

std::shared_ptr<Map> MeshSet::createMapTo(std::shared_ptr<MeshSet> to_set, std::string name, INT dim, INT* data)
{
    std::shared_ptr<Map> map(new Map(std::move(name), this->shared_from_this(), to_set, dim, data));
    this->maps.push_back(map);

    return map;    
}

//*************************************************************************************************

ParticleSet::ParticleSet(std::string name, std::shared_ptr<MeshSet> cells_set, int size)
    : Set(std::move(name), size), cells_set(cells_set)
{
    std::cout << "At constructor of [MeshSet] " << this->name << " size: " << this->size << std::endl;
}

ParticleSet::~ParticleSet()
{
    std::cout << "At distructor of [MeshSet] " << this->name << std::endl;
}

void ParticleSet::print() 
{
    std::cout << this->name << " [MeshSet] size: " << this->size << std::endl; 
}

void ParticleSet::increaseParticleCount(UINT add_count)
{
    this->add_count = add_count;
    // reserve the data structures depending on the parallelism
}

//*************************************************************************************************

Map::Map(std::string name, std::shared_ptr<Set> from_set, std::shared_ptr<Set> to_set, INT dimension, INT* map)
{
    std::cout << "At constructor of [Map] " << name << std::endl;

    if (dimension <= 0) 
    {
        printf("Map error -- negative/zero dimension for map %s\n", name.c_str());
        exit(-1);
    }
    if (map == nullptr) 
    {
        printf("Map error -- map provided is a nullptr for map %s\n", name.c_str());
        exit(-1);
    }

    this->from_set = from_set;
    this->to_set = to_set;
    this->dimension = dimension;
    this->name = std::move(name);

    this->map = util::getVecFromPtr<INT>(map, dimension * from_set->getSize());

    if (SIM->OPP_maps_base_index == 1) // convert map to 0 based indexing -- i.e. reduce each map value by 1
    {
        for (int i = 0; i < from_set->getSize() * dimension; i++)
            (this->map[i])--;
    }
}

Map::~Map()
{
    std::cout << "At distructor of [Map] " << name << std::endl;
}

Arg Map::createArg(opp_accessor acc, bool map_with_cell_index)
{
    return Arg(std::shared_ptr<BaseDat>(nullptr), 0, this->shared_from_this(), this->dimension, "int", acc, OPP_ARG_MAP, map_with_cell_index);
}

void Map::print() 
{
    std::cout << name << " [Set] size: " << dimension << std::endl; 
}

//*************************************************************************************************

Arg::Arg(std::shared_ptr<BaseDat> dat, INT idx, std::shared_ptr<Map> map, INT dim, std::string type, opp_accessor acc, opp_arg_type argtype, bool map_with_cell_index)
    : dat(dat), map(map), idx(idx), dim(dim), type(type), acc(acc), map_with_cell_index(map_with_cell_index), argtype(argtype)
{
    if (dat != nullptr)
    {
        if (type == "int")
            i_dat = &((std::dynamic_pointer_cast<Dat<INT>>(dat))->getData());
        else if (type == "double")
            r_dat = &((std::dynamic_pointer_cast<Dat<REAL>>(dat))->getData());
        else if (type == "bool")
            b_dat = &((std::dynamic_pointer_cast<Dat<BOOL>>(dat))->getData());
        else
        {
            printf("Arg error -- Data Type %s is not implemented for dat args\n", type.c_str());
            exit(-1);
        }
    }   
    if (map != nullptr)
    {
        i_map = &(map->getMapData());
    } 

    std::cout << "At constructor of [Arg]" << std::endl;
}

Arg::~Arg()
{
    std::cout << "At distructor of [Arg]" << std::endl;
}

