#include <petscksp.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "config.h"

/**
 * @brief Opaque PETSc FEM solver context.
 *
 * Heap-allocated by init_petsc_fem(). Callers store its address as a
 * @c long (Fortran: @c integer(kind=8)) so the pointer survives across
 * language boundaries without exposing internal fields.
 */
typedef struct {
    Mat      K_mat;           /**< Global stiffness matrix. */
    Vec      rhs_vec;         /**< Global right-hand side (force) vector. */
    Vec      sol_vec;         /**< Global solution (displacement) vector. */
    KSP      ksp_solver;      /**< KSP linear solver context. */

    // local solution vec 
    Vec      local_sol_vec;   /**< Local solution vector for scatter/gather operations. */
    VecScatter sol_scatter; // The MPI communication context


    PetscInt local_owned_dofs;/**< Number of DOFs owned by this MPI rank. */
    PetscInt nglob;     /**< Total global DOF count. */

    // method to solve linear system, e.g. "cg", "gmres", etc.  Set via PETSc options, e.g.:
    //   -ksp_type cg -pc_type bjacobi

    // local connectivity for one element 
    PetscInt *ibool_local;

    // specfem arrays backup 
    int num_neighbors, max_nibool;
    const int *neighbor_ranks;
    const int *npts_neighbor;
    const int *send_points; 

    // mesh info
    int nspec, NGLL3, NDIM;

    // owning rank per node, length nglob, set by build_owner_rank()
    int myrank;
    const int *owner_rank;

} PetscFemCtx;

static uint64_t 
get_hash_score(uint64_t hash_val, uint64_t rank)
{
    uint64_t score = (hash_val * UINT64_C(0x9E3779B97F4A7C15)) ^
                     (rank * UINT64_C(0xBF58476D1CE4E5B9));

    score ^= score >> 30;
    score *= UINT64_C(0xBF58476D1CE4E5B9);

    score ^= score >> 27;
    score *= UINT64_C(0x94D049BB133111EB);

    score ^= score >> 31;
    return score;
}

static double hypot3(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

/**
 * @brief Get the xyz minmax across all MPI slice
 * 
 * @param nspec Number of spectral elements.
 * @param NGLL3 Number of GLL points per element.
 * @param nglob Total number of global nodes.
 * @param ibool Connectivity array (1-based node indices, column-major).
 * @param xstore Array of x-coordinates.
 * @param ystore Array of y-coordinates.
 * @param zstore Array of z-coordinates.
 * @param minmax_out Output array for min/max values [x_min, x_max, y_min, y_max, z_min, z_max].
 * @param min_dist Output for minimum node distance
 */
static void 
get_xyz_minmax(
    int nspec, int NGLL3,
    int nglob, const int *ibool,
    const float *xstore, const float *ystore, const float *zstore,
    double *minmax_out, double *min_dist
)
{
    double minmax_loc[6] = {
        (double)xstore[0], (double)xstore[0],
        (double)ystore[0], (double)ystore[0],
        (double)zstore[0], (double)zstore[0]
    };
    for(int i = 1; i < nglob; i++) {
        if (xstore[i] < minmax_loc[0]) minmax_loc[0] = xstore[i];
        if (xstore[i] > minmax_loc[1]) minmax_loc[1] = xstore[i];
        if (ystore[i] < minmax_loc[2]) minmax_loc[2] = ystore[i];
        if (ystore[i] > minmax_loc[3]) minmax_loc[3] = ystore[i];
        if (zstore[i] < minmax_loc[4]) minmax_loc[4] = zstore[i];
        if (zstore[i] > minmax_loc[5]) minmax_loc[5] = zstore[i];
    } 

    // mpi reduce
    double gmin,gmax;
    for(int dim = 0; dim < 3; dim++) {
        MPI_Allreduce(&minmax_loc[dim*2], &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&minmax_loc[dim*2+1], &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        minmax_out[dim*2] = gmin;
        minmax_out[dim*2+1] = gmax;
    }
    
    // get NGLL on each edge 
    int NGLL = cbrt(NGLL3 + 0.001);
    if(NGLL * NGLL * NGLL != NGLL3) {
        printf("Error: NGLL3=%d is not a perfect cube, cannot determine NGLL for edge length calculation.\n", NGLL3);
        exit(1);
    }

    // loop each element to find the minimum node distance, 
    double dmin = (minmax_out[1] - minmax_out[0]) * 1.0e5; // start with a large number
    for(int e = 0; e < nspec; e++) {
        for(int iz = 0; iz < NGLL; iz ++){
        for(int iy = 0; iy < NGLL; iy ++){
        for(int ix = 0; ix < NGLL; ix ++){
            int igll3 = (iz * NGLL + iy) * NGLL + ix;
            int inode = ibool[e * NGLL3 + igll3] - 1; // 1-based to 0-based
            double x0 = xstore[inode], y0 = ystore[inode], z0 = zstore[inode];
            if(iz > 0){
                int idz = (iz-1) * NGLL * NGLL + iy * NGLL + ix;
                int inodez = ibool[e * NGLL3 + idz] - 1;
                double dz = hypot3(x0-xstore[inodez],y0-ystore[inodez],z0-zstore[inodez]);
                dmin = dmin < dz ? dmin : dz;
            }

            if(iy > 0){
                int idy = (iz * NGLL + (iy-1)) * NGLL + ix;
                int inodey = ibool[e * NGLL3 + idy] - 1;
                double dy = hypot3(x0-xstore[inodey],y0-ystore[inodey],z0-zstore[inodey]);
                dmin = dmin < dy ? dmin : dy;
            }

            if(ix > 0){
                int idx = (iz * NGLL + iy) * NGLL + (ix-1);
                int inodex = ibool[e * NGLL3 + idx] - 1;
                double dx = hypot3(x0-xstore[inodex],y0-ystore[inodex],z0-zstore[inodex]);
                dmin = dmin < dx ? dmin : dx;
            }   
        }}}
    }

    // mpi reduce to get global minimum edge length
    MPI_Allreduce(&dmin, min_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

/**
 * @brief Get the coordinate hash 
 * 
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param z Z-coordinate.
 * @param minmax_xyz Array of min/max values [x_min, x_max, y_min, y_max, z_min, z_max].
 * @param min_dist Minimum distance between nodes.
 * @return uint64_t Hash value.
 */
static uint64_t 
get_coordinate_hash(double x, double y, double z, const double *minmax_xyz, double min_dist)
{
    // quantize coordinates to a fixed precision grid based on global min/max and a predefined number of bins (e.g., 1 million)
    uint64_t num_binsx = (minmax_xyz[1] - minmax_xyz[0]) / min_dist * 10; // 10x finer than minimum edge length
    uint64_t num_binsy = (minmax_xyz[3] - minmax_xyz[2]) / min_dist * 10;
    uint64_t num_binsz = (minmax_xyz[5] - minmax_xyz[4]) / min_dist * 10;
    uint64_t x_bin = (uint64_t)((x - minmax_xyz[0]) / (minmax_xyz[1] - minmax_xyz[0]) * num_binsx);
    uint64_t y_bin = (uint64_t)((y - minmax_xyz[2]) / (minmax_xyz[3] - minmax_xyz[2]) * num_binsy);
    uint64_t z_bin = (uint64_t)((z - minmax_xyz[4]) / (minmax_xyz[5] - minmax_xyz[4]) * num_binsz);

    uint64_t hash_val = (x_bin << 42) | (y_bin << 21) | z_bin; // pack into a single 64-bit integer
    return hash_val;
}

/** @brief Cast the opaque @c long handle back to a typed pointer. */
#define CTX(handle)  ((PetscFemCtx *)(*(handle)))

/**
 * @brief Create the local-to-global mapping for the PETSc vectors/matrices.
 * 
 * @param ctx The PETSc FEM context.
 * @param l2g_out Output array for the local-to-global mapping.
 * @param rstart Output start index for this rank.
 * @param rend Output end index for this rank.
 */
static void 
create_local2global_mapping(PetscFemCtx *ctx,PetscInt *l2g_out,
                            PetscInt *rstart, PetscInt *rend) 
{
    // 1. Compute the starting Global ID for this rank
    int64_t global_offset = 0;
    int64_t local_owned = (int64_t)ctx->local_owned_dofs;
    MPI_Exscan(&local_owned, &global_offset, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (ctx->myrank == 0) global_offset = 0; // Exscan result is undefined on rank 0

    *rstart = (PetscInt)global_offset;
    *rend = (PetscInt)(global_offset + local_owned);

    // 2. Build the Local-to-Global mapping array
    int64_t *l2g_indices = (int64_t *)malloc(ctx->nglob * sizeof(int64_t));
    int64_t owned_counter = 0;
    for (int64_t i = 0; i < (int64_t) ctx->nglob; i++) {
        if (ctx->owner_rank[i] == ctx->myrank) {
            l2g_indices[i] = global_offset + owned_counter;
            owned_counter += 1;
        } 
        else {
            l2g_indices[i] = -1; // Placeholder for ghosts
        }
    }

    // sanity check if petscint is int
    if(sizeof(PetscInt) < sizeof(int64_t)) {
        for (int64_t i = 0; i < (int64_t) ctx->nglob; i++) {
            if (l2g_indices[i] > INT_MAX) {
                fprintf(stderr, "Error: Global DOF index %" PRId64 " exceeds PETSc's maximum integer limit. Consider using a 64-bit PETSc build.\n", l2g_indices[i]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // buffers
    int64_t *buf_sd = (int64_t *)malloc(ctx->num_neighbors * ctx->max_nibool * sizeof(int64_t));
    int64_t *buf_rv = (int64_t *)malloc(ctx->num_neighbors * ctx->max_nibool * sizeof(int64_t));
    MPI_Request *req_sd = (MPI_Request *)malloc(ctx->num_neighbors * sizeof(MPI_Request));
    MPI_Request *req_rv = (MPI_Request *)malloc(ctx->num_neighbors * sizeof(MPI_Request));

    for (int i = 0; i < ctx->num_neighbors; i++) {
        int npts = ctx->npts_neighbor[i];

        // Pack the Global IDs we know (the ones we own)
        for (int j = 0; j < npts; j++) {
            int local_node = ctx->send_points[i * ctx->max_nibool + j] - 1;
            buf_sd[i * ctx->max_nibool + j] = l2g_indices[local_node];
        }

        MPI_Isend(buf_sd + i * ctx->max_nibool, 
                    npts, MPI_INT64_T, ctx->neighbor_ranks[i], 99,
                    MPI_COMM_WORLD,&req_sd[i]);
        MPI_Irecv(buf_rv + i * ctx->max_nibool, 
                    npts, MPI_INT64_T, ctx->neighbor_ranks[i], 99,
                    MPI_COMM_WORLD, &req_rv[i]);
    }

    // wait finish 
    MPI_Waitall(ctx->num_neighbors, req_rv, MPI_STATUSES_IGNORE);

    // 4. Update ghost nodes with the IDs received from their actual owners
    for (int i = 0; i < ctx->num_neighbors; i++) {
        for (int j = 0; j < ctx->npts_neighbor[i]; j++) {
            int local_node = ctx->send_points[i * ctx->max_nibool + j] - 1;
            int neighbor_rank = ctx->neighbor_ranks[i];
            // Only update if we don't own it (if it's a ghost)
            if (ctx->owner_rank[local_node] == neighbor_rank) {
                l2g_indices[local_node] = buf_rv[i * ctx->max_nibool + j];
            }
        }
    }
    MPI_Waitall(ctx->num_neighbors, req_sd, MPI_STATUSES_IGNORE);

    // 5. Copy to output and clean up
    for (int i = 0; i < ctx->nglob; i++) {
        l2g_out[i] = (PetscInt)l2g_indices[i];
    }

    // cleanup 
    free(l2g_indices);
    free(buf_sd);
    free(buf_rv);
    free(req_sd);
    free(req_rv);
}

/**
 * @brief Create a node2elem adjncy list object
 * 
 * @param ctx PETSCFEMCTX context object containing mesh info
 * @param ibool connectivity array (1-based node indices, column-major)
 * @param node_elem_list array to store the node-to-element adj
 */
static void 
create_node2elem_adjncy_list(
    PetscFemCtx *ctx,
    const int *ibool,
    int **xadj_nodes,
    int **adjncy_elements
)
{
    int *num_elems_per_node = (int *)calloc(ctx->nglob, sizeof(int));
    for(int e = 0; e < ctx->nspec; e++) {
        for(int i = 0; i < ctx->NGLL3; i++) {
            int node = ibool[i + e * ctx->NGLL3] - 1;
            num_elems_per_node[node]++;
        }
    }

    // node to eleme ptr
    int *xadj = (int *)malloc((ctx->nglob + 1) * sizeof(int));
    xadj[0] = 0;
    for(int i = 0; i < ctx->nglob; i++) {
        xadj[i + 1] = xadj[i] + num_elems_per_node[i];
    }
    *xadj_nodes = xadj;

    // fill list
    int * adj_list = (int *)malloc(xadj[ctx->nglob] * sizeof(int));
    int *current_pos = (int *)calloc(ctx->nglob, sizeof(int));

    for(int e = 0; e < ctx->nspec; e++) {
        for(int i = 0; i < ctx->NGLL3; i++) {
            int node = ibool[i + e * ctx->NGLL3] - 1;
            int pos = xadj[node] + current_pos[node];
            adj_list[pos] = e;
            current_pos[node]++;
        }
    }
    *adjncy_elements = adj_list;

    free(num_elems_per_node);
    free(current_pos);
}

/**
 * @brief Build a CSR node adjacency graph mirroring test.py:create_node_adjacency.
 *
 * Adjacency uniqueness is tracked in local-node space, while stored entries
 * are the mapped PETSc global node ids from @p l2g_indices.
 */
static void
create_node_adjacency(
    PetscFemCtx *ctx,
    const int *ibool,
    int **node_adj_ptr,
    PetscInt **node_adj_data)
{
    int *node_elem_ptr = NULL;
    int *node_elem_data = NULL;
    int *node_mask = NULL;
    int *adj_ptr = NULL;
    PetscInt *adj_data = NULL;

    create_node2elem_adjncy_list(ctx, ibool, &node_elem_ptr, &node_elem_data);

    adj_ptr = (int *)calloc(ctx->nglob + 1, sizeof(int));
    node_mask = (int *)malloc(ctx->nglob * sizeof(int));
    for (int i = 0; i < ctx->nglob; i++) node_mask[i] = -1;

    for (int iglob = 0; iglob < ctx->nglob; iglob++) {
        adj_ptr[iglob + 1] = adj_ptr[iglob];
        for (int idx = node_elem_ptr[iglob]; idx < node_elem_ptr[iglob + 1]; idx++) {
            int ispec = node_elem_data[idx];
            for (int ipt = 0; ipt < ctx->NGLL3; ipt++) {
                int iglob_nb = ibool[ipt + ispec * ctx->NGLL3] - 1;
                if (node_mask[iglob_nb] != iglob) {
                    adj_ptr[iglob + 1] += 1;
                    node_mask[iglob_nb] = iglob;
                }
            }
        }
    }

    adj_data = (PetscInt *)malloc((adj_ptr[ctx->nglob] > 0 ? adj_ptr[ctx->nglob] : 1) * sizeof(PetscInt));
    for (int i = 0; i < ctx->nglob; i++) node_mask[i] = -1;

    for (int iglob = 0; iglob < ctx->nglob; iglob++) {
        int icount = 0;
        for (int idx = node_elem_ptr[iglob]; idx < node_elem_ptr[iglob + 1]; idx++) {
            int ispec = node_elem_data[idx];
            for (int ipt = 0; ipt < ctx->NGLL3; ipt++) {
                int iglob_nb = ibool[ipt + ispec * ctx->NGLL3] - 1;
                if (node_mask[iglob_nb] != iglob) {
                    adj_data[adj_ptr[iglob] + icount] = iglob_nb;
                    node_mask[iglob_nb] = iglob;
                    icount += 1;
                }
            }
        }
    }

    *node_adj_ptr = adj_ptr;
    *node_adj_data = adj_data;

    free(node_elem_ptr);
    free(node_elem_data);
    free(node_mask);
}

typedef struct {
    PetscInt *values;
    int count;
    int capacity;
} EdgeSet;

static void edge_set_init(EdgeSet *set, int initial_capacity)
{
    set->count = 0;
    set->capacity = initial_capacity > 0 ? initial_capacity : 1;
    set->values = (PetscInt *)malloc((size_t)set->capacity * sizeof(PetscInt));
}

static void edge_set_destroy(EdgeSet *set)
{
    free(set->values);
    set->values = NULL;
    set->count = 0;
    set->capacity = 0;
}

static void edge_set_add(EdgeSet *set, PetscInt value)
{
    for (int i = 0; i < set->count; i++) {
        if (set->values[i] == value) return;
    }

    if (set->count == set->capacity) {
        set->capacity *= 2;
        set->values = (PetscInt *)realloc(set->values, (size_t)set->capacity * sizeof(PetscInt));
    }

    set->values[set->count++] = value;
}

static void
create_preallocation_nnz(
    PetscFemCtx *ctx,
    const PetscInt *l2g_indices,
    PetscInt rstart,
    PetscInt rend,
    const int *node_adj_ptr,
    const PetscInt *node_adj_data,
    PetscInt **d_nnz_out,
    PetscInt **o_nnz_out)
{
    PetscInt local_owned = (PetscInt)ctx->local_owned_dofs;
    EdgeSet *owned_edges = (EdgeSet *)malloc((size_t)(local_owned > 0 ? local_owned : 1) * sizeof(EdgeSet));
    PetscInt *d_nnz = (PetscInt *)calloc((size_t)(local_owned > 0 ? local_owned : 1), sizeof(PetscInt));
    PetscInt *o_nnz = (PetscInt *)calloc((size_t)(local_owned > 0 ? local_owned : 1), sizeof(PetscInt));
    int *send_sizes = (int *)calloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1), sizeof(int));
    int *recv_sizes = (int *)calloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1), sizeof(int));
    MPI_Request *req_sd = (MPI_Request *)malloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1) * sizeof(MPI_Request));
    MPI_Request *req_rv = (MPI_Request *)malloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1) * sizeof(MPI_Request));
    PetscInt **send_bufs = NULL;
    PetscInt **recv_bufs = NULL;

    for (PetscInt row_id = 0; row_id < local_owned; row_id++) {
        edge_set_init(&owned_edges[row_id], ctx->NGLL3 * 2);
    }

    for (int iglob = 0; iglob < ctx->nglob; iglob++) {
        if (ctx->owner_rank[iglob] != ctx->myrank) continue;

        PetscInt row_id = l2g_indices[iglob] - rstart;
        for (int idx = node_adj_ptr[iglob]; idx < node_adj_ptr[iglob + 1]; idx++) {
            PetscInt global_neighbor = l2g_indices[node_adj_data[idx]];
            edge_set_add(&owned_edges[row_id], global_neighbor);
        }
    }

    for (int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        for (int ipt = 0; ipt < ctx->npts_neighbor[ineigh]; ipt++) {
            int iglob = ctx->send_points[ineigh * ctx->max_nibool + ipt] - 1;
            if (ctx->owner_rank[iglob] == neigh_rank) {
                send_sizes[ineigh] += 2 * (node_adj_ptr[iglob + 1] - node_adj_ptr[iglob]);
            }
        }

        MPI_Isend(&send_sizes[ineigh], 1, MPI_INT, neigh_rank, 100, MPI_COMM_WORLD, &req_sd[ineigh]);
        MPI_Irecv(&recv_sizes[ineigh], 1, MPI_INT, neigh_rank, 100, MPI_COMM_WORLD, &req_rv[ineigh]);
    }

    if (ctx->num_neighbors > 0) {
        MPI_Waitall(ctx->num_neighbors, req_rv, MPI_STATUSES_IGNORE);
        MPI_Waitall(ctx->num_neighbors, req_sd, MPI_STATUSES_IGNORE);
    }

    send_bufs = (PetscInt **)malloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1) * sizeof(PetscInt *));
    recv_bufs = (PetscInt **)malloc((size_t)(ctx->num_neighbors > 0 ? ctx->num_neighbors : 1) * sizeof(PetscInt *));

    for (int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        int ptr = 0;

        send_bufs[ineigh] = (PetscInt *)malloc((size_t)(send_sizes[ineigh] > 0 ? send_sizes[ineigh] : 1) * sizeof(PetscInt));
        recv_bufs[ineigh] = (PetscInt *)malloc((size_t)(recv_sizes[ineigh] > 0 ? recv_sizes[ineigh] : 1) * sizeof(PetscInt));

        for (int ipt = 0; ipt < ctx->npts_neighbor[ineigh]; ipt++) {
            int iglob = ctx->send_points[ineigh * ctx->max_nibool + ipt] - 1;
            if (ctx->owner_rank[iglob] == neigh_rank) {
                PetscInt target_global = l2g_indices[iglob];
                for (int idx = node_adj_ptr[iglob]; idx < node_adj_ptr[iglob + 1]; idx++) {
                    PetscInt neighbor_global = l2g_indices[node_adj_data[idx]];
                    send_bufs[ineigh][ptr++] = target_global;
                    send_bufs[ineigh][ptr++] = neighbor_global;
                }
            }
        }

        MPI_Isend(send_bufs[ineigh], send_sizes[ineigh], MPIU_INT, neigh_rank, 101, MPI_COMM_WORLD, &req_sd[ineigh]);
        MPI_Irecv(recv_bufs[ineigh], recv_sizes[ineigh], MPIU_INT, neigh_rank, 101, MPI_COMM_WORLD, &req_rv[ineigh]);
    }

    if (ctx->num_neighbors > 0) {
        MPI_Waitall(ctx->num_neighbors, req_rv, MPI_STATUSES_IGNORE);
        MPI_Waitall(ctx->num_neighbors, req_sd, MPI_STATUSES_IGNORE);
    }

    for (int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        for (int idx = 0; idx < recv_sizes[ineigh]; idx += 2) {
            PetscInt target_global = recv_bufs[ineigh][idx];
            PetscInt neighbor_global = recv_bufs[ineigh][idx + 1];
            PetscInt row_id = target_global - rstart;
            edge_set_add(&owned_edges[row_id], neighbor_global);
        }
    }

    for (PetscInt row_id = 0; row_id < local_owned; row_id++) {
        for (int idx = 0; idx < owned_edges[row_id].count; idx++) {
            PetscInt global_neighbor = owned_edges[row_id].values[idx];
            if (global_neighbor >= rstart && global_neighbor < rend) d_nnz[row_id] += 1;
            else o_nnz[row_id] += 1;
        }
    }

    for (PetscInt row_id = 0; row_id < local_owned; row_id++) {
        edge_set_destroy(&owned_edges[row_id]);
    }
    free(owned_edges);
    for (int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        free(send_bufs[ineigh]);
        free(recv_bufs[ineigh]);
    }
    free(send_bufs);
    free(recv_bufs);
    free(send_sizes);
    free(recv_sizes);
    free(req_sd);
    free(req_rv);

    *d_nnz_out = d_nnz;
    *o_nnz_out = o_nnz;
}


/**
 * @brief Create and preallocate the global stiffness matrix and vectors.
 *
 * Uses @c ctx->local_owned_dofs (set by build_owner_rank()) as the local
 * row/column size.  The caller must supply the diagonal and off-diagonal
 * non-zero counts so PETSc can allocate memory efficiently without
 * dynamic reallocation during assembly.
 *
 * @param[in] ctx           Solver context handle.
 * @param [in] ibool         Element connectivity (1-based node indices, column-major).
 */
static void 
setup_petsc_impl(PetscFemCtx *ctx, const int *ibool)
{
    PetscInt local_owned = (PetscInt)ctx->local_owned_dofs;
    PetscInt rstart, rend;
    PetscInt *l2g_indices = (PetscInt *)malloc(ctx->nglob * sizeof(PetscInt));

    // 1. Create Local-to-Global mapping and get ownership range
    create_local2global_mapping(ctx, l2g_indices, &rstart, &rend);

    // =========================================================================
    // 2. Build local node adjacency, mirroring test.py:create_node_adjacency
    // =========================================================================
    int *node_adj_ptr = NULL;
    PetscInt *node_adj_data = NULL;
    create_node_adjacency(ctx, ibool, &node_adj_ptr, &node_adj_data);

    PetscInt *d_nnz = NULL;
    PetscInt *o_nnz = NULL;
    create_preallocation_nnz(ctx, l2g_indices, rstart, rend, node_adj_ptr, node_adj_data, &d_nnz, &o_nnz);

    // =========================================================================
    // 6. Create PETSc Objects
    // =========================================================================
    PetscInt local_dofs = local_owned * ctx->NDIM;

    PetscCallVoid(MatCreate(PETSC_COMM_WORLD, &ctx->K_mat));
    PetscCallVoid(MatSetType(ctx->K_mat, MATMPIBAIJ));
    PetscCallVoid(MatSetSizes(ctx->K_mat, local_dofs, local_dofs, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCallVoid(MatSetBlockSize(ctx->K_mat, ctx->NDIM));
    PetscCallVoid(MatMPIBAIJSetPreallocation(ctx->K_mat, ctx->NDIM, 0, d_nnz, 0, o_nnz));

    ISLocalToGlobalMapping l2g_mapping;
    PetscCallVoid(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, ctx->NDIM, ctx->nglob, 
                                               l2g_indices, PETSC_COPY_VALUES, &l2g_mapping));
    PetscCallVoid(MatSetLocalToGlobalMapping(ctx->K_mat, l2g_mapping, l2g_mapping));
    
    // 
    PetscCallVoid(VecCreateMPI(PETSC_COMM_WORLD, local_dofs, PETSC_DETERMINE, &ctx->rhs_vec));
    PetscCallVoid(VecSetBlockSize(ctx->rhs_vec, ctx->NDIM));
    PetscCallVoid(VecSetLocalToGlobalMapping(ctx->rhs_vec, l2g_mapping));
    PetscCallVoid(VecDuplicate(ctx->rhs_vec, &ctx->sol_vec));

    // local solution vector
    PetscCallVoid(VecCreateSeq(PETSC_COMM_SELF, ctx->nglob * ctx->NDIM, &ctx->local_sol_vec));
    PetscCallVoid(VecSetBlockSize(ctx->local_sol_vec, ctx->NDIM));

    // scatter info
    IS from_is, to_is;
    PetscCallVoid(ISCreateBlock(PETSC_COMM_SELF, ctx->NDIM, ctx->nglob, 
                                l2g_indices, PETSC_COPY_VALUES, &from_is));
    PetscInt *seq_indices = (PetscInt *)malloc(ctx->nglob * sizeof(PetscInt));
    for(int i = 0; i < ctx->nglob; i++) {
        seq_indices[i] = i;
    }
    PetscCallVoid(ISCreateBlock(PETSC_COMM_SELF, ctx->NDIM, ctx->nglob, 
                                seq_indices, PETSC_COPY_VALUES, &to_is));
    PetscCallVoid(VecScatterCreate(ctx->sol_vec, from_is, ctx->local_sol_vec, to_is, &ctx->sol_scatter));

    // 5. Cleanup the temporary Index Sets
    PetscCallVoid(ISDestroy(&from_is));
    PetscCallVoid(ISDestroy(&to_is));
    PetscCallVoid(ISLocalToGlobalMappingDestroy(&l2g_mapping));

    // Cleanup
    free(l2g_indices); free(d_nnz); free(o_nnz);
    free(node_adj_ptr); free(node_adj_data);
    free(seq_indices);
}

static void 
setup_petsc_impl1(PetscFemCtx *ctx, const int *ibool)
{
    PetscInt local_owned = (PetscInt)ctx->local_owned_dofs;
    PetscInt rstart, rend;
    PetscInt *l2g_indices = (PetscInt *)malloc(ctx->nglob * sizeof(PetscInt));

    // 1. Create Local-to-Global mapping and get ownership range
    create_local2global_mapping(ctx, l2g_indices, &rstart, &rend);

    // =========================================================================
    // 2. Build local node adjacency, mirroring test.py:create_node_adjacency
    // =========================================================================
    int *node_adj_ptr = NULL;
    PetscInt *node_adj_data = NULL;
    create_node_adjacency(ctx, ibool, &node_adj_ptr, &node_adj_data);

    PetscInt *d_nnz = NULL;
    PetscInt *o_nnz = NULL;
    create_preallocation_nnz(ctx, l2g_indices, rstart, rend, node_adj_ptr, node_adj_data, &d_nnz, &o_nnz);

    // =========================================================================
    // 6. Create PETSc Objects
    // =========================================================================
    PetscInt local_dofs = local_owned * ctx->NDIM;

    PetscCallVoid(MatCreate(PETSC_COMM_WORLD, &ctx->K_mat));
    PetscCallVoid(MatSetType(ctx->K_mat, MATMPIBAIJ));
    PetscCallVoid(MatSetSizes(ctx->K_mat, local_dofs, local_dofs, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCallVoid(MatSetBlockSize(ctx->K_mat, ctx->NDIM));
    PetscCallVoid(MatMPIBAIJSetPreallocation(ctx->K_mat, ctx->NDIM, 0, d_nnz, 0, o_nnz));
    
    // ** CRITICAL FIX ** Tells PETSc the incoming Fortran dense element arrays are Column-Major!
    //PetscCallVoid(MatSetOption(ctx->K_mat, MAT_ROW_ORIENTED, PETSC_FALSE)); 

    ISLocalToGlobalMapping l2g_mapping;
    PetscCallVoid(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, ctx->NDIM, ctx->nglob, 
                                               l2g_indices, PETSC_COPY_VALUES, &l2g_mapping));
    PetscCallVoid(MatSetLocalToGlobalMapping(ctx->K_mat, l2g_mapping, l2g_mapping));

    PetscCallVoid(VecCreateMPI(PETSC_COMM_WORLD, local_dofs, PETSC_DETERMINE, &ctx->rhs_vec));
    PetscCallVoid(VecSetBlockSize(ctx->rhs_vec, ctx->NDIM));
    PetscCallVoid(VecSetLocalToGlobalMapping(ctx->rhs_vec, l2g_mapping));
    PetscCallVoid(VecDuplicate(ctx->rhs_vec, &ctx->sol_vec));
    PetscCallVoid(ISLocalToGlobalMappingDestroy(&l2g_mapping));

    // Cleanup
    free(l2g_indices); free(d_nnz); free(o_nnz);
    free(node_adj_ptr); free(node_adj_data);
}

/**
 * @brief Determine DOF ownership using the "Lowest Rank Wins" rule.
 *
 * For every node shared across MPI partitions the rank with the smallest
 * rank ID is declared the owner.  After the call, @p owner_rank[i] holds
 * the owning rank for node @c i, and @c ctx->local_owned_dofs is set to
 * the number of nodes owned by @p myrank.
 *
 * @p ibool_interfaces_ext_mesh is stored in column-major (Fortran) order
 * with leading dimension @p max_nibool:
 * @code
 *   node = ibool_interfaces_ext_mesh[ipoin + max_nibool * iinterface]  // 1-based
 * @endcode
 *
 * @param[in]  h                            Solver context handle.
 * @param[in]  nglob                        Total nodes on this partition.
 * @param[in]  myrank                       MPI rank of this process.
 * @param[in]  num_interfaces_ext_mesh      Number of shared interfaces.
 * @param[in]  nibool_interfaces_ext_mesh   Node count per interface (length @p num_interfaces_ext_mesh).
 * @param[in]  ibool_interfaces_ext_mesh    1-based node indices, column-major layout.
 * @param[in]  max_nibool                   Leading dimension of @p ibool_interfaces_ext_mesh.
 * @param[in]  my_neighbors_ext_mesh        MPI rank of each interface neighbour (length @p num_interfaces_ext_mesh).
 * @param[in]  xstore,ystore,zstore         Node coordinates (length @p nglob).
 * @param[in]  ibool                        Element connectivity (1-based node indices, column-major
 * @param[out] owner_rank                   Owning rank per node (length @p nglob).
 */
void FC_FUNC_(setup_petsc,SETUP_PETSC)(long *h,
                                    int nglob, int myrank,
                                    int nspec,int NGLL3, int NDIM,
                                    int num_interfaces_ext_mesh,
                                    const int *nibool_interfaces_ext_mesh,
                                    const int *ibool_interfaces_ext_mesh,
                                    int  max_nibool,
                                    const int *my_neighbors_ext_mesh,
                                    const float *xstore,
                                    const float *ystore,
                                    const float *zstore,
                                    const int *ibool,
                                    int  *owner_rank)
{
    PetscFemCtx *ctx = (PetscFemCtx *)malloc(sizeof(PetscFemCtx));
    *h = (long)ctx;
    int i, iinterface, ipoin, inode, neighbor;

    if(myrank == 0) {
        printf("-------------------------------\n");
        printf("creating ownerships ...\n");
        printf("-------------------------------\n");
    }

    /* Initialise: this rank claims every node it touches. */
    for (i = 0; i < nglob; i++)
        owner_rank[i] = myrank;

    // create global_x/y/z min/max
    double minmax_xyz[6],min_dist;
    get_xyz_minmax(nspec,NGLL3,nglob,ibool, xstore, ystore, zstore, minmax_xyz, &min_dist);

    // copy a reference of ghost node information to the context for later use in setup_petsc() when we need to compute the non-zero pattern of the matrix
    ctx->num_neighbors = num_interfaces_ext_mesh;
    ctx->max_nibool = max_nibool;
    ctx->neighbor_ranks = my_neighbors_ext_mesh;
    ctx->npts_neighbor = nibool_interfaces_ext_mesh;
    ctx->send_points = ibool_interfaces_ext_mesh;

    /* Knockout loop: a neighbour with a lower rank steals ownership. */
    for(int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        int num_points = ctx->npts_neighbor[ineigh];
        const int *points = &ctx->send_points[ineigh * ctx->max_nibool];

        for(int ipoin = 0; ipoin < num_points; ipoin++) {
            int node = points[ipoin] - 1; // convert to 0-based index

            // get coordinate hash for this node
            double x = (double)xstore[node];
            double y = (double)ystore[node];
            double z = (double)zstore[node];
            uint64_t hash_val = get_coordinate_hash(x, y, z, minmax_xyz, min_dist);
            uint64_t my_score = get_hash_score(hash_val, (uint64_t)owner_rank[node]);
            uint64_t neigh_score = get_hash_score(hash_val, (uint64_t)neigh_rank);

            if(neigh_score > my_score) {
                owner_rank[node] = neigh_rank;
            }
        }
    }

    ctx->local_owned_dofs = 0;
    for (i = 0; i < nglob; i++)
        if (owner_rank[i] == myrank)
            ctx->local_owned_dofs++;

    printf("Rank %d claiming ownership of all %d nodes of %d total nodes\n", myrank, (int)ctx->local_owned_dofs, nglob);

    ctx->nglob = nglob;

    // backup of owner_rank in the context for later use in setup_petsc() when we need to compute the non-zero pattern of the matrix
    ctx->owner_rank = owner_rank;
    ctx->myrank = myrank;
    ctx->nspec = nspec;
    ctx->NGLL3 = NGLL3;
    ctx->NDIM = NDIM;

    // allocate local connectivity for one element, which will be used in setup_petsc() to compute the non-zero pattern of the matrix
    ctx->ibool_local = (PetscInt *)malloc(NGLL3 * sizeof(PetscInt));

    // petsc interface
    if(1 == 1) {
        setup_petsc_impl(ctx, ibool);
    }
    else {
        setup_petsc_impl1(ctx, ibool);
    }
}

/**
 * @brief Accumulate one element's stiffness block into the global matrix.
 *
 * @param[in] h              Solver context handle.
 * @param[in] local_indices 0-based global DOF indices (length @p dofs_per_elem).
 * @param[in] elem_matrix    Dense element stiffness values (row-major).
 */
void FC_FUNC_(fill_mat_petsc,FILL_MAT_PETSC)(long *h,
                                             const int *local_indices,
                                             const PetscScalar *k_elem)
{
    PetscFemCtx   *ctx = CTX(h);
    PetscErrorCode  ierr;
    PetscInt i;
    int NGLL3 = ctx->NGLL3;
    int NDIM = ctx->NDIM;

    // check if petscint is int
    PetscInt *id_input = NULL;
    if(sizeof(PetscInt) > sizeof(int)) {
        id_input = ctx->ibool_local;
        for(i = 0; i < (PetscInt) NGLL3; i++) {
            int node_id = local_indices[i];
            id_input[i] = (PetscInt) node_id;
        }
    }
    else{
        id_input = (PetscInt *)local_indices;
    }

    // 3. Insert the entire 24x24 block using just the 8 Node IDs!
    PetscCallVoid(
        MatSetValuesBlockedLocal(ctx->K_mat, 
                                NGLL3, id_input, // Row nodes
                                NGLL3, id_input, // Column nodes
                                k_elem,                // Pointer to dense matrix data
                                ADD_VALUES)
    );

}

/**
 * @brief Accumulate one element's force vector into the global RHS vector.
 * @param[in] h              Solver context handle.
 * @param[in] local_indices 0-based global DOF indices (length @p NGLL3 * @p NDIM).
 * @param[in] f_elem         Dense element force values (length @p NGLL3 * @p NDIM, ordered by DOF then node, i.e. [u1_node1, v1_node1, w1_node1, u
 */
void FC_FUNC_(fill_vec_petsc,FILL_VEC_PETSC)(long *h,
                                             const int *local_indices,
                                             const PetscScalar *f_elem)
{
    PetscFemCtx   *ctx = CTX(h);
    PetscErrorCode  ierr;
    int NGLL3 = ctx->NGLL3;

    // 3. Insert the entire NGLL3*NDIM-value array using just the NGLL3 Node IDs!
    // ADD_VALUES is critical here so shared node forces accumulate!
    PetscInt *id_input = NULL;
    if(sizeof(PetscInt) > sizeof(int)) {
        id_input = ctx->ibool_local;
        for(int i = 0; i < NGLL3; i++) {
            int node_id = local_indices[i];
            id_input[i] = (PetscInt) node_id;
        }
    }
    else{
        id_input = (PetscInt *)local_indices;
    }

    PetscCallVoid(
        VecSetValuesBlockedLocal(ctx->rhs_vec, 
                               NGLL3, id_input, // The NGLL3 node IDs
                               f_elem,            // The NGLL3*NDIM force values
                               ADD_VALUES)
    );
}

/**
 * @brief Finalise global matrix and RHS vector assembly.
 *
 * Triggers the underlying MPI communication that moves off-rank contributions
 * to their owning processes.  Must be called after all assemble_element()
 * calls and before solve_system().
 *
 * @param[in] h  Solver context handle.
 */
void FC_FUNC_(assemble_petsc,ASSEMBLE_PETSC)(long *h)
{
    PetscFemCtx   *ctx = CTX(h);

    PetscCallVoid(MatAssemblyBegin(ctx->K_mat, MAT_FINAL_ASSEMBLY));
    PetscCallVoid(MatAssemblyEnd  (ctx->K_mat, MAT_FINAL_ASSEMBLY));

    PetscCallVoid(VecAssemblyBegin(ctx->rhs_vec));
    PetscCallVoid(VecAssemblyEnd  (ctx->rhs_vec));

}

/**
 * @brief Solve the assembled linear system  \f$ K \, u = f \f$.
 *
 * Creates a KSP solver, attaches the stiffness matrix as the system operator,
 * then calls KSPSolve().  The solver type and preconditioner can be
 * overridden at run-time via PETSc command-line options, e.g.:
 * @code
 *   -ksp_type cg -pc_type bjacobi
 * @endcode
 *
 * @param[in] h  Solver context handle.
 */
void FC_FUNC_(solve_petsc,SOLVE_PETSC)(long *h)
{
    PetscFemCtx   *ctx = CTX(h);

    PetscCallVoid(KSPCreate(PETSC_COMM_WORLD, &ctx->ksp_solver));
    PetscCallVoid(KSPSetOperators(ctx->ksp_solver, ctx->K_mat, ctx->K_mat));
    PetscCallVoid(KSPSetFromOptions(ctx->ksp_solver));
    PetscCallVoid(KSPSolve(ctx->ksp_solver, ctx->rhs_vec, ctx->sol_vec));
}

/**
 * @brief Copy owned solution DOFs from PETSc into the caller's native array and broadcast to ghosts.
 *
 * Only entries for which @p owner_array[i] is non-zero are extracted from PETSc.
 * After extraction, this routine performs a non-blocking MPI exchange to update 
 * all ghost nodes across the partition boundaries.
 *
 * @param[in]     h                   Solver context handle.
 * @param[in,out] displ               Output native array (size ctx->nglob * NDIM).
 */
void FC_FUNC_(extract_petsc,EXTRACT_PETSC)(long *h, double *displ)
{
    PetscFemCtx       *ctx = CTX(h);
    int NDIM = ctx->NDIM;

#ifndef USE_PETSC_SCATTER_SOL_VEC 

    // 1. Fire the MPI communication! 
    // This fetches all ghost node values from neighbor ranks instantly.
    PetscCallVoid(VecScatterBegin(ctx->sol_scatter, ctx->sol_vec, ctx->local_sol_vec, INSERT_VALUES, SCATTER_FORWARD));
    PetscCallVoid(VecScatterEnd  (ctx->sol_scatter, ctx->sol_vec, ctx->local_sol_vec, INSERT_VALUES, SCATTER_FORWARD));

    // 2. Extract the raw array and copy to Fortran
    const PetscScalar *array;
    PetscCallVoid(VecGetArrayRead(ctx->local_sol_vec, &array));

    // Copy to the flat Fortran array
    for (int i = 0; i < ctx->nglob * NDIM; i++) {
        displ[i] = (double)array[i];
    }

    PetscCallVoid(VecRestoreArrayRead(ctx->local_sol_vec, &array));

#else 

    // =========================================================================
    // 1. Extract owned data from PETSc
    // =========================================================================
    const PetscScalar *petsc_data;
    PetscCallVoid(VecGetArrayRead(ctx->sol_vec, &petsc_data));

    int petsc_node_idx = 0; // Tracker for our position in the packed PETSc array
    for (int i = 0; i < ctx->nglob; i++) {
        // We only have data for nodes we strictly own
        if (ctx->owner_rank[i] == ctx->myrank) {
            // Unpack the NDIM block (e.g., X, Y, Z)
            for (int dim = 0; dim < NDIM; dim++) {
                displ[i * NDIM + dim] = petsc_data[petsc_node_idx * NDIM + dim];
            }
            // Advance our tracker in the packed PETSc array
            petsc_node_idx++;
        }
    }
    PetscCallVoid(VecRestoreArrayRead(ctx->sol_vec, &petsc_data));

    // =========================================================================
    // 2. Allocate buffers for MPI communication
    // =========================================================================
    size_t buf_size = (size_t)ctx->num_neighbors * ctx->max_nibool * NDIM * sizeof(double);
    double *send_buffer = (double *)malloc(buf_size);
    double *recv_buffer = (double *)malloc(buf_size);  
    MPI_Request *req_sd = (MPI_Request *)malloc((size_t)ctx->num_neighbors * sizeof(MPI_Request));
    MPI_Request *req_rv = (MPI_Request *)malloc((size_t)ctx->num_neighbors * sizeof(MPI_Request));

    // =========================================================================
    // 3. Pack data and Initiate Non-Blocking Sends/Receives
    // =========================================================================
    for(int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        int num_points = ctx->npts_neighbor[ineigh];
        const int *points = &ctx->send_points[ineigh * ctx->max_nibool];

        // Packing Loop
        for(int ipt = 0; ipt < num_points; ipt++) {
            int node = points[ipt] - 1; // convert to 0-based index
            bool is_owned = (ctx->owner_rank[node] == ctx->myrank);
            
            for(int idim = 0; idim < NDIM; idim++) {
                size_t offset = (size_t)ineigh * ctx->max_nibool * NDIM + ipt * NDIM + idim;
                
                if (is_owned) {
                    send_buffer[offset] = displ[node * NDIM + idim];
                } else {
                    send_buffer[offset] = 0.0; // Dummy value, neighbor will ignore
                }
            }
        }

        // FIRE MPI COMMUNICATORS OUTSIDE THE POINT LOOP!
        // We must pass the starting memory address of this neighbor's specific block
        size_t start_offset = (size_t)ineigh * ctx->max_nibool * NDIM;
        
        MPI_Isend(&send_buffer[start_offset], num_points * NDIM, MPI_DOUBLE, 
                  neigh_rank, 0, MPI_COMM_WORLD, &req_sd[ineigh]);
                  
        MPI_Irecv(&recv_buffer[start_offset], num_points * NDIM, MPI_DOUBLE, 
                  neigh_rank, 0, MPI_COMM_WORLD, &req_rv[ineigh]);
    }

    // =========================================================================
    // 4. Wait for receives to finish and unpack
    // =========================================================================
    MPI_Waitall(ctx->num_neighbors, req_rv, MPI_STATUSES_IGNORE);

    for(int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        int num_points = ctx->npts_neighbor[ineigh];
        const int *points = &ctx->send_points[ineigh * ctx->max_nibool];
        
        for(int ipt = 0; ipt < num_points; ipt++) {
            int node = points[ipt] - 1; 
            
            for(int idim = 0; idim < NDIM; idim++) {
                size_t offset = (size_t)ineigh * ctx->max_nibool * NDIM + ipt * NDIM + idim;
                
                // Only accept data if it came from the true owner!
                if (ctx->owner_rank[node] == neigh_rank) {
                    displ[node * NDIM + idim] = recv_buffer[offset];
                }
            }
        }
    }

    // =========================================================================
    // 5. Clean up
    // =========================================================================
    // Wait for sends to finish before freeing the send_buffer memory
    MPI_Waitall(ctx->num_neighbors, req_sd, MPI_STATUSES_IGNORE);

    free(send_buffer);
    free(recv_buffer);
    free(req_sd);
    free(req_rv);
#endif
}

/**
 * @brief Destroy all PETSc objects, finalise the library, and free the context.
 *
 * Sets @c *h to @c 0 after freeing so the caller's handle is nulled out and
 * any accidental second call is detectable.
 *
 * @param[in,out] h  Solver context handle; set to @c 0 on return.
 */
void FC_FUNC_(cleanup_petsc,CLEANUP_PETSC)(long *h)
{
    PetscFemCtx   *ctx = CTX(h);

    PetscCallVoid(KSPDestroy(&ctx->ksp_solver));
    PetscCallVoid(VecDestroy(&ctx->rhs_vec));
    PetscCallVoid(VecDestroy(&ctx->sol_vec));
    PetscCallVoid(MatDestroy(&ctx->K_mat));
    free(ctx->ibool_local);
    PetscCallVoid(VecDestroy(&ctx->local_sol_vec));
    PetscCallVoid(VecScatterDestroy(&ctx->sol_scatter));

    free(ctx);
    *h = 0;
}
