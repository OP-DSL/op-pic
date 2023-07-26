
#pragma once

#include <oppic_lib.h>
#include <list>
#include <cmath>
#include "opp_defs.h"

namespace opp {

    class CellApproximator {
    
    public:
        
        CellApproximator(const opp_dat nodePos_dat, double gridSpacing); 

        virtual ~CellApproximator();

        // CellApproximator(const CellApproximator& obj) = delete;
        // static CellApproximator* getInstance();
        // static void releaseInstance();

        opp_move_status getCellIndexKernel(const double *point_pos, int* current_cell_index, double* point_lc,  
            const double *cell_volume, const double *cell_det, const int *cell_connectivity);
        
        const std::array<opp_point, 2>& generateBoundingBox(const opp_dat node_pos_dat);

        const std::array<opp_point, 2>& getBoundingBox();

        opp_point getCentroidOfBox(const opp_point& coordinate);

        const std::vector<int>& generateStructMeshToCellIndexVec(const opp_dat cell_volume_dat, const opp_dat cell_det_dat, 
            const opp_map cell_connectivity_map);

        const std::vector<int>& getStructMeshToCellIndexVec();

        bool isCoordinateInBoundingBox(const opp_point& point);

        int findClosestCellIndex(const opp_point& targetPosition);
        
        const std::vector<opp_point>& generateStructCoordinateVec();

        const std::vector<opp_point>& getStructuredCoordinateVec();

        const opp_ipoint& calculateGridDimensions();
        
        const opp_ipoint& getGridDimensions();

        const std::vector<unsigned int>& getHopCountsVec();
        
        void countHopsFromVec(const std::string& name, const std::vector<size_t>& vec);

        void move(const opp_dat pos_dat, opp_dat cell_index_dat, opp_dat lc_dat, const opp_dat cell_volume_dat, 
            const opp_dat cell_det_dat, const opp_map cell_connectivity_map);

        size_t getCellIndexMappingIndex(const opp_point& position);
    private:

        // CellApproximator* instancePtr = nullptr;

        double gridSpacing = 0.0;
        double oneOverGridSpacing = 0.0;
        std::array<opp_point,2> boundingBox; // index 0 is min, index 1 is max
        opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
        opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);
        opp_ipoint gridDimensions;
        std::vector<int> structMeshToCellIndexMap;
        std::vector<opp_point> coordinateVec;
        std::vector<unsigned int> hopCountsVec;

        const double ONE_OVER_SIX = (1.0 / 6.0);
        const int N_PER_C = 4;
        const int DET_FIELDS = 4;
        const int NEIGHB_C = 4;
    };
};

