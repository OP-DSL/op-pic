#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

#define INT int
#define UINT unsigned int
#define REAL double

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
    
    
    private:
        std::string name, loves;

        Simulation();   
        virtual ~Simulation();

        static Simulation* instancePtr;
    };

};

    

