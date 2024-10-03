#pragma once

#include <limits>
#include <array>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include <cstring>
#include <inttypes.h>
#include <stdio.h>
#include <stdarg.h>
#include <memory>
#include <cmath>
#include <sys/queue.h>
#include <sys/time.h>
#include <random>
#include <unistd.h>

#ifdef USE_MPI
    namespace opp {
        class Comm;
        class GlobalParticleMover;
    };
#else
    #define Comm int
    #define GlobalParticleMover int
#endif

#ifdef USE_PETSC
    #include <petscksp.h>
#endif

#ifdef USE_THRUST
    #include <thrust/device_vector.h>
    #include <thrust/host_vector.h>
    #define THRUST_REAL thrust::device_vector<double>
    #define THRUST_INT thrust::device_vector<int>
#elif defined(USE_SYCL)
    #include <sycl/sycl.hpp>
    #define THRUST_REAL void
    #define THRUST_INT void    
    extern sycl::queue* opp_queue;
#else
    #define THRUST_REAL void
    #define THRUST_INT void
#endif

#include <opp_params.h>
#include <opp_profiler.h>

constexpr double MAX_REAL = std::numeric_limits<double>::max();
constexpr double MIN_REAL = std::numeric_limits<double>::min();

constexpr int MAX_INT = std::numeric_limits<int>::max();
constexpr int MIN_INT = std::numeric_limits<int>::min();

#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

#define OPP_RUN_ON_ROOT(command) if (OPP_rank == OPP_ROOT command)

#define ASSIGN_CENTROID_TO_DIM(K)                                   \
    if (coordinate.K + this->gridSpacing <= maxCoordinate.K) {      \
        centroid.K = coordinate.K + this->gridSpacing * 0.5;         \
    }                                                               \
    else {                                                          \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5;        \
    }                                                               \

#define CHECK(cmd) { int err = cmd; if (err != MPI_SUCCESS) opp_abort(std::to_string(err)); }

//*************************************************************************************************
#ifdef DEBUG_LOG
    #define OPP_DBG true
#else
    #define OPP_DBG false
#endif

#ifdef OPP_ROOT
    #undef OPP_ROOT
#endif
#define OPP_ROOT 0

#define MAX_CELL_INDEX     INT_MAX

#define UNUSED(expr) do { (void)(expr); } while (0)

#define OPP_READ        0
#define OPP_WRITE       1
#define OPP_RW          2
#define OPP_INC         3
#define OPP_MIN         4
#define OPP_MAX         5

#define OPP_ARG_GBL     0
#define OPP_ARG_DAT     1
#define OPP_ARG_MAP     2

#define ZERO_double    0.0
#define ZERO_float     0.0f
#define ZERO_int       0
#define ZERO_uint      0
#define ZERO_ll        0
#define ZERO_ull       0
#define ZERO_bool      0
#define OPP_REAL_ZERO  0.0
#define OPP_INT_ZERO   0
#define OPP_BOOL_ZERO  0

#define OPP_DEFAULT_GPU_THREADS_PER_BLOCK 32

#define LOG_STR_LEN 100000

//*************************************************************************************************
enum opp_iterate_type
{
    OPP_ITERATE_ALL = 1,
    OPP_ITERATE_INJECTED,
};

enum opp_move_status 
{
    OPP_MOVE_DONE = 0,
    OPP_NEED_MOVE,
    OPP_NEED_REMOVE,
};

enum opp_data_type 
{
    DT_INT = 0,
    DT_REAL,
};

enum DeviceType
{
    Device_CPU = 1,
    Device_GPU = 2,
};

enum Dirty
{
    NotDirty = 0,
    Device = 1,
    Host = 2,
};

enum opp_mapping
{
    OPP_Map_Default = 0,
    OPP_Map_from_Mesh_Rel,
    OPP_Map_from_Inj_part,
    OPP_Map_to_Mesh_Rel,
};

enum opp_reset
{
    OPP_Reset_Core = 0,
    OPP_Reset_Set,
    OPP_Reset_ieh,
    OPP_Reset_inh,
    OPP_Reset_All,
};

enum opp_reduc_comm
{
    OPP_Reduc_NO_Comm = 0,
    OPP_Reduc_SUM_Comm,
    OPP_Reduc_MIN_Comm,
    OPP_Reduc_MAX_Comm,
};

struct part_index {
    int start;
    int end;
};

enum opp_fill_type {
    OPP_HoleFill_All = 0,
    OPP_Sort_All,
    OPP_Shuffle_All,
    OPP_Sort_Periodic,
    OPP_Shuffle_Periodic,
};

//*************************************************************************************************
struct opp_set_core;
typedef struct opp_set_core *opp_set;

struct opp_map_core;
typedef struct opp_map_core *opp_map;

struct opp_dat_core;
typedef struct opp_dat_core *opp_dat;

typedef int opp_access;       /* holds OPP_READ, OPP_WRITE, OPP_RW, OPP_INC, OPP_MIN, OPP_MAX */
typedef int opp_arg_type;     /* holds OPP_ARG_GBL, OPP_ARG_DAT, OPP_ARG_MAP */

struct opp_arg {
    int index;                  /* index */
    opp_dat dat;                /* dataset */
    opp_map map;                /* indirect mapping */
    opp_dat p2c_map;            /* double indirect mapping - used for only particles */
    int dim;                    /* dimension of data */
    int idx;                    /* size (for sequential execution) */
    int size;                   /* size (for sequential execution) */
    char *data;                 /* data on host */
    char *data_d;               /* data on device (for CUDA execution) */
    int *map_data;              /* data on host */
    int *map_data_d;            /* data on device (for CUDA execution) */
    char const *type;           /* datatype */
    opp_access acc;             /* opp_accessor OPP_READ, OPP_WRITE, OPP_RW, OPP_INC, OPP_MIN, OPP_MAX */
    opp_arg_type argtype;
    int sent;                   /* flag to indicate if this argument has data in flight under non-blocking MPI comms*/
    int opt;                    /* flag to indicate if this argument is in use */
    opp_mapping mesh_mapping;
};

struct opp_set_core {
    int index;                              /* index */
    int size;                               /* number of elements in set */
    char const *name;                       /* name of set */
    int core_size;                          /* number of core elements in an mpi process*/
    int exec_size;                          /* number of additional imported elements to be executed */
    int nonexec_size;                       /* number of additional imported elements that are not executed */

    bool is_particle;                       /* is it a particle set */
    int set_capacity;                       /* capacity of the allocated array */
    int diff;                               /* number of particles to change */
    int particle_size;                      /* size of particle */
    std::vector<int>* indexes_to_remove;
    opp_dat mesh_relation_dat = NULL;
    std::vector<opp_dat>* particle_dats;
    std::map<int, part_index>* cell_index_v_part_index_map;
    // int* particle_statuses;
    // int* particle_statuses_d;
    int particle_remove_count;
    int* particle_remove_count_d;
    void* mpi_part_buffers;
    opp_set cells_set;
};

struct opp_map_core {
    int index;                  /* index */
    opp_set from;             /* set pointed from */
    opp_set to;               /* set pointed to */
    int dim;                    /* dimension of pointer */
    int *map;                   /* array defining pointer */
    int *map_d;                 /* device array defining pointer */
    char const *name;           /* name of pointer */
    int user_managed;           /* indicates whether the user is managing memory */
    opp_dat p2c_dat;
};

struct opp_dat_core {
    int index;                  /* index */
    opp_set set;                /* set on which data is defined */
    int dim;                    /* dimension of data */
    int size;                   /* size of each element in dataset */
    char *data;                 /* data on host */
    char *data_d;               /* data on device (GPU) */
    char *data_swap_d;          /* data on device (GPU) - used for swapping */
    char const *type;           /* datatype */
    char const *name;           /* name of dataset */
    char *buffer_d;             /* buffer for MPI halo sends on the devidce */
    char *buffer_d_r;           /* buffer for MPI halo receives on the devidce */
    int dirtybit;               /* flag to indicate MPI halo exchange is needed*/
    Dirty dirty_hd;             /* flag to indicate dirty status on host and device */
    int user_managed;           /* indicates whether the user is managing memory */
    void *mpi_buffer;           /* ponter to hold the mpi buffer struct for the opp_dat*/    
    void *mpi_reduc_buffer;     /* ponter to hold the mpi reduction buffer struct for the opp_dat*/ 
    opp_reduc_comm reduc_comm;  /* flag to check whether the dat is in between reduction communication */

    std::vector<char*>* thread_data;
    bool is_cell_index;

    THRUST_INT *thrust_int;
    THRUST_REAL *thrust_real;

    THRUST_INT *thrust_int_sort;
    THRUST_REAL *thrust_real_sort;

    opp_map p2c_map;
};

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

constexpr double opp_zero_double16[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
constexpr int opp_zero_int16[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//*******************************************************************************
namespace opp {

    //*******************************************************************************
    class BoundingBox {

    public:
        BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate);
        BoundingBox(const opp_dat node_pos_dat, int dim);
        ~BoundingBox();

        const opp_point& getLocalMin() const;
        const opp_point& getLocalMax() const;
        const opp_point& getGlobalMin() const;
        const opp_point& getGlobalMax() const;
        bool isCoordinateInBoundingBox(const opp_point& point);
        bool isCoordinateInGlobalBoundingBox(const opp_point& point);
        inline int getDim() const { return dim; }

    private:
        void generateGlobalBoundingBox(int count);
        int dim = 0;

        std::array<opp_point,2> boundingBox; // index 0 is min, index 1 is max
        std::array<opp_point,2> globalBoundingBox; // index 0 is min, index 1 is max
    };

    //*******************************************************************************
    class GlobalToLocalCellIndexMapper {
    
    public:
        //*******************************************************************************
        // This will contain mappings for halo indices too
        GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat);
        virtual ~GlobalToLocalCellIndexMapper();

        int map(const int globalIndex);
        
    private:
        std::map<int,int> globalToLocalCellIndexMap;
    };

    //*******************************************************************************
    class CellMapper {
    
    public:
        CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, 
            const std::shared_ptr<Comm> comm = nullptr);
        ~CellMapper();

        opp_point getCentroidOfBox(const opp_point& coordinate);
        size_t findStructuredCellIndex2D(const opp_point& position);
        size_t findStructuredCellIndex3D(const opp_point& position);
        int findClosestCellIndex(const size_t& structCellIdx);
        int findClosestCellRank(const size_t& structCellIdx);
        void reduceInterNodeMappings(int callID);
        void convertToLocalMappings(const opp_dat global_cell_id_dat);
        void enrichStructuredMesh(const int index, const int cell_index, const int rank = 0);
        void printStructuredMesh(const std::string msg, int *array, size_t size, bool printToFile = true);
        void createStructMeshMappingArrays();
        void waitBarrier();
        void lockWindows();
        void unlockWindows();

    // private:
        const std::shared_ptr<BoundingBox> boundingBox = nullptr;
        const double gridSpacing = 0.0;
        const double oneOverGridSpacing = 0.0;
        const opp_point& minGlbCoordinate;  
        const std::shared_ptr<Comm> comm = nullptr;

        opp_ipoint globalGridDims;
        opp_ipoint localGridStart, localGridEnd;

        size_t globalGridSize = 0;
        
        int* structMeshToCellMapping = nullptr;         // This contain mapping to local cell indices
        int* structMeshToRankMapping = nullptr;         // This contain mapping to residing mpi rank

#ifdef USE_MPI        
        MPI_Win win_structMeshToCellMapping;
        MPI_Win win_structMeshToRankMapping;
#endif
    };
};

extern std::shared_ptr<opp::BoundingBox> boundingBox;
extern std::shared_ptr<opp::CellMapper> cellMapper;
extern bool useGlobalMove;

#ifdef USE_MPI 
    extern std::shared_ptr<opp::Comm> comm;
    extern std::unique_ptr<opp::GlobalParticleMover> globalMover;
#else
    extern std::shared_ptr<Comm> comm;
    extern std::unique_ptr<GlobalParticleMover> globalMover;
#endif

extern int OPP_hybrid_gpu;
extern double OPP_hybrid_balance;
extern int OPP_maps_base_index;
extern int OPP_auto_soa;
extern int OPP_gpu_direct;
extern double OPP_part_alloc_mult;
// extern int OPP_auto_sort;
extern int OPP_fill_period;
extern opp_fill_type OPP_fill_type;
extern int OPP_mpi_part_alloc_mult;
extern int OPP_rank;
extern int OPP_comm_size;
extern int OPP_comm_iteration;
extern int OPP_max_comm_iteration;
extern int OPP_iter_start;
extern int OPP_iter_end;
extern int *OPP_mesh_relation_data;
extern int *OPP_mesh_relation_data_d;
extern int OPP_main_loop_iter;
extern int OPP_gpu_threads_per_block;
extern size_t OPP_gpu_shared_mem_per_block;
extern int OPP_part_cells_set_size;
extern int OPP_part_comm_count_per_iter;
extern int OPP_move_max_hops;

extern std::vector<opp_set> opp_sets;
extern std::vector<opp_map> opp_maps;
extern std::vector<opp_dat> opp_dats;

extern std::unique_ptr<opp::Params> opp_params;
extern std::unique_ptr<opp::Profiler> opp_profiler;
