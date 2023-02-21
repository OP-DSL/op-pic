#include "opp_simulation.h"

using namespace opp;

//*************************************************************************************************

Simulation* Simulation::instancePtr = NULL;

Simulation::Simulation()
{
  std::cout << "At constructor of [Simulation] " << std::endl;
}

Simulation::~Simulation()
{
  std::cout << "At distructor of [Simulation] " << std::endl;
}

Simulation* Simulation::getInstance()
{
  if (instancePtr == NULL)
  {   
    instancePtr = new Simulation();

    std::cout << "Creating new instance of [Simulation] " << std::endl;
    return instancePtr;
  }

  std::cout << "Returning old instance of [Simulation] " << std::endl;
  return instancePtr;
}

void Simulation::releaseInstance()
{
  if (instancePtr != NULL)
  {   
    delete instancePtr;
    instancePtr = NULL;

    std::cout << "released instance of [Simulation] " << std::endl;
  }
}

std::shared_ptr<MeshSet> Simulation::createMeshSet(std::string name, int size)
{
  return nullptr;
}

std::shared_ptr<ParticleSet> Simulation::createParticleSet(std::string name, std::shared_ptr<MeshSet> cells_set, int size)
{
  return nullptr;
}

UINT Simulation::registerDat(std::shared_ptr<BaseDat> dat)
{
  opp_dats.push_back(dat);
  return (opp_dats.size() - 1);
}

UINT Simulation::registerMap(std::shared_ptr<Map> map)
{
  opp_maps.push_back(map);
  return (opp_maps.size() - 1);
}

UINT Simulation::registerSet(std::shared_ptr<Set> set)
{
  opp_sets.push_back(set);
  return (opp_sets.size() - 1);
}

void Simulation::setValues(std::string name, std::string loves)
{
  this->name = name;
  this->loves = loves;
}

void Simulation::print()
{
  std::cout << name << " Loves " <<
          loves << "." << std::endl;
}

//*************************************************************************************************
