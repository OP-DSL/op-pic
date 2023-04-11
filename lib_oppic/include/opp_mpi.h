#pragma once

#include <mpi.h>
#include <oppic_lib.h>
#include <cmath>
#include <sys/queue.h>

#ifdef OPP_MPI_ROOT
    #undef OPP_MPI_ROOT
#endif
#define OPP_MPI_ROOT 0





/** extern variables for halo creation and exchange**/
extern MPI_Comm OP_MPI_WORLD;
extern MPI_Comm OP_MPI_GLOBAL;

/*******************************************************************************
* MPI halo list data type
*******************************************************************************/
typedef struct 
{
    op_set set;  // set related to this list
    int size; // number of elements in this list 
    int *ranks; // MPI ranks to be exported to or imported from 
    int ranks_size; // number of MPI neighbors to be exported to or imported from
    int *disps; // displacements for the starting point of each rank's element list 
    int *sizes; // number of elements exported to or imported from each ranks 
    int *list; // the list of all elements
} halo_list_core;

typedef halo_list_core *halo_list;

/*******************************************************************************
* Data structures related to MPI level partitioning
*******************************************************************************/
typedef struct 
{  
    op_set set; // set to which this partition info blongs to  
    int *g_index;  // global index of each element held in this MPI process 
    int *elem_part;  // partition to which each element belongs  
    int is_partitioned; // indicates if this set is partitioned 1 if partitioned 0 if not
} part_core;

typedef part_core *part;

/*******************************************************************************
* Data structure to hold mpi communications of an op_dat
*******************************************************************************/
#define NAMESIZE 20
typedef struct 
{ 
    char name[NAMESIZE]; // name of this op_dat 
    int size; // size of this op_dat 
    int index; // index of this op_dat
    int count; // total number of times this op_dat was halo exported
    int bytes; // total number of bytes halo exported for this op_dat in this kernel
} op_dat_mpi_comm_info_core;

typedef op_dat_mpi_comm_info_core *op_dat_mpi_comm_info;

/*******************************************************************************
* Buffer struct used in non-blocking mpi halo sends/receives
*******************************************************************************/
typedef struct 
{ 
    char *buf_exec; // buffer holding exec halo to be exported; 
    char *buf_nonexec; // buffer holding nonexec halo to be exported; 
    MPI_Request *s_req; // pointed to hold the MPI_Reqest for sends  
    MPI_Request *r_req; // pointed to hold the MPI_Reqest for receives  
    int s_num_req; // number of send MPI_Reqests in flight at a given time for this op_dat 
    int r_num_req; // number of receive MPI_Reqests in flight at a given time for this op_dat
} op_mpi_buffer_core;

typedef op_mpi_buffer_core *op_mpi_buffer;

/*******************************************************************************
* Buffer struct used in non-blocking mpi particle sends/receives
*******************************************************************************/
struct opp_mpi_part_buffer 
{ 
    char *buf_import;
    int buf_import_capacity;
    int buf_import_index;       // not used
    char *buf_export;
    int buf_export_capacity;    // init to -1
    int buf_export_index;       // init to zero
};

struct opp_all_mpi_part_buffers 
{ 
    int total_recv;
    std::map<int,opp_mpi_part_buffer> buffers;
    std::map<int,int> import_counts; // rank -> count
    std::map<int,int> export_counts; // rank -> count
    std::vector<MPI_Request> send_req;
    std::vector<MPI_Request> recv_req;
    std::vector<int> neighbours;
};

struct opp_particle_comm_data
{ 
    int cell_residing_rank;
    int local_index;

};

/** external variables **/
extern int OP_part_index;
extern part *OP_part_list;
extern int **orig_part_range;

/** export list on the device **/
extern int **export_exec_list_d;
extern int **export_exec_list_disps_d;
extern int **export_nonexec_list_d;
extern int **export_nonexec_list_disps_d;
extern int **export_nonexec_list_partial_d;
extern int **import_nonexec_list_partial_d;
extern int *set_import_buffer_size;
extern int **import_exec_list_disps_d;
extern int **import_nonexec_list_disps_d;

extern halo_list *OP_export_exec_list; // EEH list
extern halo_list *OP_import_exec_list; // IEH list

extern halo_list *OP_import_nonexec_list; // INH list
extern halo_list *OP_export_nonexec_list; // ENH list

extern int OP_part_index;
extern part *OP_part_list;
extern int **orig_part_range;

extern std::map<int, std::map<int, opp_particle_comm_data>> opp_part_comm_neighbour_data;

/*******************************************************************************
* Data Type to hold sliding planes info
*******************************************************************************/

typedef struct 
{
    int index;
    int coupling_group_size;
    int *coupling_proclist;

    int num_ifaces;
    int *iface_list;

    int *nprocs_per_int;
    int **proclist_per_int;
    int **nodelist_send_size;
    int ***nodelist_send;

    int max_data_size;
    char ***send_buf;
    MPI_Request **requests;
    MPI_Status **statuses;

    char *OP_global_buffer;
    int OP_global_buffer_size;

    int gbl_num_ifaces;
    int *gbl_iface_list;
    int *nprocs_per_gint;
    int **proclist_per_gint;

    int gbl_offset;
    op_map cellsToNodes;
    op_dat coords;
    op_dat mark;
} op_export_core;

typedef op_export_core *op_export_handle;

typedef struct 
{
    int index;
    int nprocs;
    int *proclist;
    int gbl_offset;
    op_dat coords;
    op_dat mark;
    int max_dat_size;
    int num_my_ifaces;
    int *iface_list;
    int *nprocs_per_int;
    int **proclist_per_int;
    int *node_size_per_int;
    int **nodelist_per_int;
    char ***recv_buf;
    int *recv2int;
    int *recv2proc;
    MPI_Request *requests;
    MPI_Status *statuses;
    double *interp_dist;
} op_import_core;

typedef op_import_core *op_import_handle;


/*******************************************************************************
* Functions declared in opp_mpi_partition.cpp
*******************************************************************************/
#if __cplusplus
extern "C" {
#endif
void opp_partition_kway(op_map primary_map);
void opp_partition_destroy();

/* static */ void partition_all(op_set primary_set, int my_rank, int comm_size);
/* static */ void renumber_maps(int my_rank, int comm_size);
/* static */ void migrate_all(int my_rank, int comm_size);

/* static */ int frequencyof(int value, int *array, int size);
/* static */ int find_mode(int *array, int size);
/* static */ int compare_all_sets(op_set target_set, op_set other_sets[], int size);
/* static */ int *create_exp_list_2(op_set set, int *temp_list, halo_list h_list, int *part_list, int size, int comm_size, int my_rank);
/* static */ void create_imp_list_2(op_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank);
/* static */ int partition_from_set(op_map map, int my_rank, int comm_size, int **part_range);
/* static */ int partition_to_set(op_map map, int my_rank, int comm_size, int **part_range);
#if __cplusplus
}
#endif

void opp_sanitize_all_maps();
void opp_desanitize_all_maps();
/*******************************************************************************/

/*******************************************************************************
* Functions declared in opp_mpi_halo.cpp
*******************************************************************************/

void opp_halo_create();
void opp_halo_destroy();

void decl_partition(op_set set, int *g_index, int *partition);
void get_part_range(int **part_range, int my_rank, int comm_size, MPI_Comm Comm);
int get_partition(int global_index, int *part_range, int *local_index, int comm_size);
int get_global_index(int local_index, int partition, int *part_range, int comm_size);
int get_global_index(int local_index, int partition, int *part_range, int comm_size);
void find_neighbors_set(halo_list List, int *neighbors, int *sizes, int *ranks_size, int my_rank, int comm_size, MPI_Comm Comm);
void create_list(int *list, int *ranks, int *disps, int *sizes, int *ranks_size, int *total, int *temp_list, 
        int size, int comm_size, int my_rank);
void create_export_list(op_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank);
void create_import_list(op_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, 
        int ranks_size, int comm_size, int my_rank);
/* static */  void create_nonexec_import_list(op_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank);
/* static */  void create_nonexec_export_list(op_set set, int *temp_list, halo_list h_list, int total_size,
        int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank);
/*******************************************************************************/

/*******************************************************************************
* Inline Utility Functions
*******************************************************************************/

inline int opp_get_uniform_local_size(int global_size, int rank) 
{
    int local_size = global_size / OPP_comm_size;
    int remainder = (int)fmod(global_size, OPP_comm_size);

    if (rank < remainder) 
    {
        local_size = local_size + 1;
    }

    return local_size;
}

inline int opp_get_uniform_local_size(int global_size) 
{
    return opp_get_uniform_local_size(global_size, OPP_my_rank);
}

template <typename T> 
inline void opp_uniform_scatter_array(T *g_array, T *l_array, int g_size, int l_size, int elem_size) 
{
    int *sendcnts = (int *)malloc(OPP_comm_size * sizeof(int));
    int *displs = (int *)malloc(OPP_comm_size * sizeof(int));
    int disp = 0;

    for (int i = 0; i < OPP_comm_size; i++) 
    {
        sendcnts[i] = elem_size * opp_get_uniform_local_size(g_size, i) * sizeof(T);
        //printf("RANK %d %d\t| sendcount %d | %d\n", OPP_my_rank, g_size, sendcnts[i], opp_get_uniform_local_size(g_size, i));
    }
    for (int i = 0; i < OPP_comm_size; i++) 
    {
        displs[i] = disp;
        disp = disp + sendcnts[i];
        //printf("RANK %d %d\t| displs %d\n", OPP_my_rank, g_size, displs[i]);
    }

    MPI_Scatterv((char*)g_array, sendcnts, displs, MPI_CHAR, 
        (char*)l_array, (l_size * elem_size * sizeof(T)), MPI_CHAR, OPP_MPI_ROOT, MPI_COMM_WORLD);

    free(sendcnts);
    free(displs);
}

/*******************************************************************************/

void opp_partition(op_set prime_set, op_map prime_map, op_dat data = NULL);
bool opp_check_all_done(oppic_set set);
void opp_wait_all_particles(oppic_set set);
bool opp_check_part_need_comm(int map0idx, oppic_set set, int particle_index);
void opp_exchange_particles(oppic_set set);
void opp_partition_core(op_set prime_set, op_map prime_map, op_dat data);
void opp_particle_comm_init();
void opp_particle_set_comm_init(oppic_set set);