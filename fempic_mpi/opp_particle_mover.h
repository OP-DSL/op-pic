#pragma once

#include <oppic_lib.h>
#include "opp_cell_mapper.h"
#include "opp_particle_mover_kernel.h"

#ifdef USE_OMP
#include <omp.h>
#endif        

namespace opp {

    class ParticleMover {
    
    public:

        //*******************************************************************************
        ParticleMover(const double gridSpacing, int dim, const opp_dat node_pos_dat, const opp_dat cellVolume_dat, 
            const opp_dat cellDet_dat, const opp_map cellConnectivity_map)
            : gridSpacing(gridSpacing), dim(dim), cellVolume_dat(cellVolume_dat), cellDet_dat(cellDet_dat), 
              cellConnectivity_map(cellConnectivity_map) {

#ifdef ENABLE_MPI
            comm = std::make_shared<Comm>(MPI_COMM_WORLD);
#endif
            boundingBox = std::make_shared<BoundingBox>(node_pos_dat, dim, comm);
            cellMapper = std::make_shared<CellMapper>(boundingBox, gridSpacing, comm);

            cellMapper->generateStructMeshToCellIndexMapping(cellVolume_dat, cellDet_dat, cellConnectivity_map);
        }

        //*******************************************************************************
        ~ParticleMover() { };

// #define USE_OMP
#ifdef USE_OMP

        //*******************************************************************************
        inline void move(opp_set set, const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat) { 
            
            opp_profiler->start("MoveApprox");

            const double* pos = (const double*)pos_dat->data;
            int* cellIndex = (int*)cellIndex_dat->data;
            double* lc = (double*)lc_dat->data;
            int set_size = set->size;
            
            // hopCountsVec.resize(size);
            // std::fill(hopCountsVec.begin(), hopCountsVec.end(), 0); 

            opp_init_particle_move(set, 0, nullptr);

            int nthreads = omp_get_max_threads();

            #pragma omp parallel for
            for (int thr = 0; thr < nthreads; thr++)
            {
                size_t start  = ((size_t)set_size * thr) / nthreads;
                size_t finish = ((size_t)set_size * (thr+1)) / nthreads;
                
                int cellIdx = MAX_CELL_INDEX;

                for (size_t i = start; i < finish; i++) {   
                    
                    opp_point* point = (opp_point*)&(pos[i * 3]);
                    cellIndex[i] = cellMapper->findClosestCellIndex(*point);

                    if (cellIndex[i] == MAX_CELL_INDEX) {
                        part_remove_count_per_thr[thr] += 1;
                        // std::cout << "Error... " << i << " " << cellIndex[i] << " " << cellVolume_dat->set->size << std::endl;
                        continue;
                    }

                    opp_move_var m;

                    do {
                        cellIdx = cellIndex[i];

                        m.move_status = getCellIndexKernel(
                            &((const double*)pos)[i * 3], 
                            &((int*)cellIndex)[i],
                            &((double*)lc)[i * 4],
                            &((double*)cellVolume_dat->data)[cellIdx], 
                            &((double*)cellDet_dat->data)[cellIdx * cellDet_dat->dim], 
                            &((int*)cellConnectivity_map->map)[cellIdx * cellConnectivity_map->dim]);
                        
                        // hopCountsVec[i]++;

                    } while (opp_part_check_status(m, cellIdx, set, i, thr, thr));
                }
            }

            opp_finalize_particle_move(set);

            opp_profiler->end("MoveApprox");
        }

        #else

        //*******************************************************************************
        inline void move(opp_set set, const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat) { 
            
            opp_profiler->start("MoveApprox");

            const double* pos = (const double*)pos_dat->data;
            int* cellIndex = (int*)cellIndex_dat->data;
            double* lc = (double*)lc_dat->data;
            int size = set->size;
            int cellIdx = MAX_CELL_INDEX;

            // hopCountsVec.resize(size);
            // std::fill(hopCountsVec.begin(), hopCountsVec.end(), 0); 

            opp_init_particle_move(set, 0, nullptr);

            for (int i = 0; i < size; i++) {   
                
                opp_point* point = (opp_point*)&(pos[i * 3]);
                cellIndex[i] = cellMapper->findClosestCellIndex(*point);

                if (cellIndex[i] == MAX_CELL_INDEX) {
                    set->particle_remove_count++;
                    // std::cout << "Error... " << i << " " << cellIndex[i] << " " << cellVolume_dat->set->size << std::endl;
                    continue;
                }

                opp_move_var m;

                do {
                    cellIdx = cellIndex[i];

                    m.move_status = getCellIndexKernel(
                        &((const double*)pos)[i * 3], 
                        &((int*)cellIndex)[i],
                        &((double*)lc)[i * 4],
                        &((double*)cellVolume_dat->data)[cellIdx], 
                        &((double*)cellDet_dat->data)[cellIdx * cellDet_dat->dim], 
                        &((int*)cellConnectivity_map->map)[cellIdx * cellConnectivity_map->dim]);
                    
                    // hopCountsVec[i]++;

                } while (opp_part_check_status(m, cellIdx, set, i, set->particle_remove_count));
            }

            opp_finalize_particle_move(set);

            opp_profiler->end("MoveApprox");
        }

        #endif

    private:

        std::shared_ptr<BoundingBox> boundingBox;
        std::shared_ptr<CellMapper> cellMapper;
        std::shared_ptr<Comm> comm;

        const opp_dat cellDet_dat;
        const opp_dat cellVolume_dat;
        const opp_map cellConnectivity_map;

        const double gridSpacing;
        const int dim;
    };
};