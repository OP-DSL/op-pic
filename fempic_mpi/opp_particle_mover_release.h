#pragma once

#include <oppic_lib.h>
#include "opp_cell_mapper.h"
#include "opp_particle_mover_kernel.h"

#ifdef ENABLE_MPI
#include "opp_global_move.h"
#endif  

const bool X_DEBUG = false; // OP_DEBUG;

#ifndef ENABLE_MPI
#define ENABLE_MPI
#endif

namespace opp {

    class ParticleMover {
    
    public:

        //*******************************************************************************
        ParticleMover(const double gridSpacing, int dim, const opp_dat node_pos_dat, const opp_dat cellVolume_dat, 
            const opp_dat cellDet_dat, const opp_dat global_cell_id_dat, const opp_map cellConnectivity_map)
            : gridSpacing(gridSpacing), dim(dim), cellVolume_dat(cellVolume_dat), cellDet_dat(cellDet_dat), 
              global_cell_id_dat(global_cell_id_dat), cellConnectivity_map(cellConnectivity_map) {
            
            opp_profiler->start("SetupMover");

            useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

            if (useGlobalMove) {
#ifdef ENABLE_MPI
                comm = std::make_shared<Comm>(MPI_COMM_WORLD);
                globalMover = std::make_unique<GlobalParticleMover>(comm->comm_parent);
#endif
                boundingBox = std::make_shared<BoundingBox>(node_pos_dat, dim, comm);

                cellMapper = std::make_shared<CellMapper>(boundingBox, gridSpacing, comm);

                cellMapper->generateStructMeshToGlobalCellMappings(cellVolume_dat, cellDet_dat, global_cell_id_dat, 
                                                                    cellConnectivity_map);

                cellMapper->generateGlobalToLocalCellIndexMapping(global_cell_id_dat);
            }

            opp_profiler->reg("GlbToLocal");
            opp_profiler->reg("GblMv_Move");
            opp_profiler->reg("GblMv_AllMv");
            for (int i = 0; i < 10; i++) {
                std::string profName = std::string("Mv_AllMv") + std::to_string(i);
                opp_profiler->reg(profName);
            }

            opp_profiler->end("SetupMover");
        }

        //*******************************************************************************
        ~ParticleMover() { };

// #define USE_OMP
#ifdef USE_OMP

        //*******************************************************************************
        inline void move(opp_set set, const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat, opp_dat part_id) { 
            
            opp_profiler->start("Move");

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

            opp_profiler->end("Move");
        }

#elif defined ENABLE_MPI

        //*******************************************************************************
        inline void move(opp_set set, const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat, opp_dat partID_dat) { 
            
            opp_profiler->start("Move");

            this->partID = partID_dat;
            this->partCellID = cellIndex_dat;

            // lambda function for multi hop particle movement
            auto multihop_mover = [&](const int i) {

                int& cellIdx = ((int*)cellIndex_dat->data)[i];

                if (cellIdx == MAX_CELL_INDEX) {
                    return;
                }

                opp_move_var m;

                do {
                    m.move_status = getCellIndexKernel(
                        &((const double*)pos_dat->data)[i * 3], 
                        &((int*)cellIndex_dat->data)[i],
                        &((double*)lc_dat->data)[i * 4],
                        &((double*)cellVolume_dat->data)[cellIdx], 
                        &((double*)cellDet_dat->data)[cellIdx * 16],        // 16 -> cellDet_dat->dim
                        &((int*)cellConnectivity_map->map)[cellIdx * 4]);   // 4 -> cellConnectivity_map->dim

                } while (opp_part_check_status(m, cellIdx, set, i, set->particle_remove_count));
            };

            // ----------------------------------------------------------------------------
            opp_init_particle_move(set, 0, nullptr);
            
            if (useGlobalMove) {
                
                opp_profiler->start("GblMv_Move");

                globalMover->initGlobalMove();

                // check whether particles needs to be moved over global move routine
                for (int i = OPP_iter_start; i < OPP_iter_end; i++) {   
                    
                    int* cellIdx = &((int*)cellIndex_dat->data)[i];
                    const opp_point* point = (const opp_point*)&(((double*)pos_dat->data)[i * 3]);

                    // check for global move, and if satisfy global move criteria, then remove the particle from current rank
                    if (checkForGlobalMove(set, *point, i, *cellIdx)) {
                        
                        *cellIdx = MAX_CELL_INDEX;
                        set->particle_remove_count++;
                        
                        continue;  
                    }
                }

                opp_profiler->end("GblMv_Move");

                globalMover->communicate(set);
            }

            opp_profiler->start("Mv_AllMv0");

            // ----------------------------------------------------------------------------
            // check whether all particles not marked for global comm is within cell, 
            // and if not mark to move between cells within the MPI rank, mark for neighbour comm
            for (int i = OPP_iter_start; i < OPP_iter_end; i++) { 
               
                multihop_mover(i);
            }

            opp_profiler->end("Mv_AllMv0");

            // ----------------------------------------------------------------------------
            // finalize the global move routine and iterate over newly added particles and check whether they need neighbour comm
            if (useGlobalMove && globalMover->finalize(set) > 0) {
                
                opp_profiler->start("GblMv_AllMv");

                { // could remove this block if the local cell index is directly stamped in cell index field
                    
                    opp_profiler->start("GlbToLocal");

                    // Change cell index from global to local. Do only for the globally moved particles 
                    for (int i = (set->size - set->diff); i < set->size; i++) { 

                        int globalCellIndex = ((int*)cellIndex_dat->data)[i];   
                        ((int*)cellIndex_dat->data)[i] = cellMapper->getLocalCellIndexFromGlobal(globalCellIndex);

                        if (((int*)cellIndex_dat->data)[i] == MAX_CELL_INDEX) {
                            opp_printf("ParticleMover", "Error... particle %d does not have a correct cell index after global move", i);
                            set->particle_remove_count++;
                        }
                    }

                    opp_profiler->end("GlbToLocal");
                }

                // check whether the new particle is within cell, and if not move between cells within the MPI rank, 
                // mark for neighbour comm. Do only for the globally moved particles 
                for (int i = (set->size - set->diff); i < set->size; i++) { 
                      
                    multihop_mover(i);                 
                }

                opp_profiler->end("GblMv_AllMv");
            }

            // ----------------------------------------------------------------------------
            // Do neighbour communication and if atleast one particle is received by the currect rank, 
            // then iterate over the newly added particles
            while (opp_finalize_particle_move(set)) {
                
                std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
                opp_profiler->start(profName);
                
                opp_init_particle_move(set, 0, nullptr);

                // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
                for (int i = OPP_iter_start; i < OPP_iter_end; i++) { 
                    
                    multihop_mover(i);
                }

                opp_profiler->end(profName);
            }

            opp_profiler->end("Move");
        }

#else
        //*******************************************************************************
        inline void move(opp_set set, const opp_dat pos_dat, opp_dat cellIndex_dat, opp_dat lc_dat, opp_dat part_id) { 
            
            opp_profiler->start("Move");

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
                
                // Since SEQ use global indices, we can simply use findClosestGlobalCellIndex
                size_t structCellIdx = cellMapper->findStructuredCellIndex(*point);
                cellIndex[i] = (int)cellMapper->findClosestGlobalCellIndex(structCellIdx); 
                
                if (cellIndex[i] == MAX_CELL_INDEX) { // Particle is outside the mesh, need to remove

                    set->particle_remove_count++;
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

            opp_profiler->end("Move");
        }

#endif

        //*******************************************************************************
        inline std::string getGlobalCellIndex(const int localCellId) {
            return std::to_string(((int*)global_cell_id_dat->data)[localCellId]);
        }

        //*******************************************************************************
        // returns true, if the current particle needs to be removed from the rank
        inline bool checkForGlobalMove(opp_set set, const opp_point& point, const int partIndex, int& cellIdx) {

#ifdef ENABLE_MPI            
            size_t structCellIdx = cellMapper->findStructuredCellIndex(point);

            if (structCellIdx == MAX_CELL_INDEX) { // This happens when point is out of the unstructured mesh
                if (OP_DEBUG)
                    opp_printf("GlobalMove", 
                    "Remove %d [Struct cell index invalid - strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                        partIndex, structCellIdx, point.x, point.y, point.z);
                return true;
            }

            int structCellRank = cellMapper->findClosestCellRank(structCellIdx);

            // Check whether the paticles need global moving, if yes start global moving process, 
            // if no, move to the closest local cell
            if (structCellRank != OPP_rank) {

                if (structCellRank == MAX_CELL_INDEX) {
                    if (OP_DEBUG)
                        opp_printf("GlobalMove", 
                        "Remove %d [Rank invalid - strCellRank:%d gblCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                            partIndex, structCellRank, cellMapper->findClosestGlobalCellIndex(structCellIdx), structCellIdx, 
                            point.x, point.y, point.z);
                    return true;
                }

                // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
                size_t globalCellIndex = cellMapper->findClosestGlobalCellIndex(structCellIdx);

                if (globalCellIndex == MAX_CELL_INDEX) {
                    if (OP_DEBUG)
                        opp_printf("GlobalMove", 
                        "Remove %d [CellIdx invalid - strCellRank:%d gblCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                            partIndex, structCellRank, globalCellIndex, structCellIdx, point.x, point.y, point.z);
                    return true;
                }

                // if the new rank is not the current rank, mark the particle to be sent via global comm
                globalMover->markParticleToMove(set, partIndex, structCellRank, globalCellIndex);

                // if (OP_DEBUG)
                //     opp_printf("GlobalMove", "Mark part %d [Move to rank %d gblCellIdx %d]", 
                //         partIndex, structCellRank, globalCellIndex);
               
                return true;
            }
            else {
                
                // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
                cellIdx = (int)cellMapper->findClosestLocalCellIndex(structCellIdx);
            }
#endif           
            return false;
        }

    private:

        std::shared_ptr<BoundingBox> boundingBox;
        std::shared_ptr<CellMapper> cellMapper;
        std::shared_ptr<Comm> comm;
        std::unique_ptr<GlobalParticleMover> globalMover;
        
        opp_dat partID = nullptr;
        opp_dat partCellID = nullptr;
        std::stringstream markedForMove;
        std::stringstream markedForRemove;

        const double gridSpacing;
        const int dim = 3;
        const opp_dat cellVolume_dat;
        const opp_dat cellDet_dat;
        const opp_dat global_cell_id_dat;
        const opp_map cellConnectivity_map;

        bool useGlobalMove = true;
    };
};