#pragma once

#include <limits>
#include <array>
#include <string>
#include <algorithm>
#include <sstream>

#ifdef ENABLE_MPI
    #include "opp_mpi.h"
#else
    #define Comm int
    #define MPI_Win int
    #define GlobalParticleMover int
#endif

#ifdef USE_OMP
#include <omp.h>
#endif        

constexpr double MAX_REAL = std::numeric_limits<double>::max();
constexpr double MIN_REAL = std::numeric_limits<double>::min();

constexpr int MAX_INT = std::numeric_limits<int>::max();
constexpr int MIN_INT = std::numeric_limits<int>::min();

struct opp_point {
    opp_point(double _x, double _y, double _z) {
        x = _x; 
        y = _y;
        z = _z;
    };
    opp_point() { };

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct opp_ipoint {
    opp_ipoint(int _x, int _y, int _z) {
        x = _x; 
        y = _y;
        z = _z;
    };
    opp_ipoint() { };

    int x = 0;
    int y = 0;
    int z = 0;
};

// void opp_abort(const std::string& s);

#define CHECK(cmd) { int err = cmd; if (err != MPI_SUCCESS) opp_abort(std::to_string(err)); }
