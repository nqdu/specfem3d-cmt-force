#include <petscksp.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>

#include "config.h"

/**
 * @file petsc_routines_coo.c
 * @brief COO-API variant of petsc_routines.c (matrix only).
 *
 * The global stiffness matrix is assembled through PETSc's coordinate (COO)
 * interface (MatSetPreallocationCOO / MatSetValuesCOO) instead of
 * MatSetValuesBlockedLocal. The RHS and solution vectors keep the original
 * VecSetValuesBlockedLocal / VecAssemblyBegin/End pipeline. The external
 * Fortran-callable API is identical to petsc_routines.c, so this file is a
 * drop-in replacement.
 */

typedef struct {
    Mat        K_mat;
    Vec        rhs_vec;
    Vec        sol_vec;
    KSP        ksp_solver;

    Vec        local_sol_vec;
    VecScatter sol_scatter;

    PetscInt   local_owned_dofs;
    PetscInt   nglob;

    PetscInt * conn; // connectivity matrix

    int        num_neighbors, max_nibool;
    const int *neighbor_ranks;
    const int *npts_neighbor;
    const int *send_points;

    int        nspec, NGLL3, NDIM;
    int        myrank;
    const int *owner_rank;
} PetscFemCtx;

#define CTX(handle)  ((PetscFemCtx *)(*(handle)))

/* =========================================================================
 * Helpers shared with petsc_routines.c — kept verbatim so this file is
 * self-contained and can be built in isolation.
 * ========================================================================= */

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

static void
get_xyz_minmax(
    int nspec, int NGLL3,
    int nglob, const int *ibool,
    const float *xstore, const float *ystore, const float *zstore,
    double *minmax_out, double *min_dist)
{
    double minmax_loc[6] = {
        (double)xstore[0], (double)xstore[0],
        (double)ystore[0], (double)ystore[0],
        (double)zstore[0], (double)zstore[0]
    };
    for (int i = 1; i < nglob; i++) {
        if (xstore[i] < minmax_loc[0]) minmax_loc[0] = xstore[i];
        if (xstore[i] > minmax_loc[1]) minmax_loc[1] = xstore[i];
        if (ystore[i] < minmax_loc[2]) minmax_loc[2] = ystore[i];
        if (ystore[i] > minmax_loc[3]) minmax_loc[3] = ystore[i];
        if (zstore[i] < minmax_loc[4]) minmax_loc[4] = zstore[i];
        if (zstore[i] > minmax_loc[5]) minmax_loc[5] = zstore[i];
    }

    double gmin, gmax;
    for (int dim = 0; dim < 3; dim++) {
        MPI_Allreduce(&minmax_loc[dim*2],   &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&minmax_loc[dim*2+1], &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        minmax_out[dim*2]     = gmin;
        minmax_out[dim*2 + 1] = gmax;
    }

    int NGLL = (int)cbrt(NGLL3 + 0.001);
    if (NGLL * NGLL * NGLL != NGLL3) {
        printf("Error: NGLL3=%d is not a perfect cube\n", NGLL3);
        exit(1);
    }

    double dmin = (minmax_out[1] - minmax_out[0]) * 1.0e5;
    for (int e = 0; e < nspec; e++) {
        for (int iz = 0; iz < NGLL; iz++) {
        for (int iy = 0; iy < NGLL; iy++) {
        for (int ix = 0; ix < NGLL; ix++) {
            int igll3 = (iz * NGLL + iy) * NGLL + ix;
            int inode = ibool[e * NGLL3 + igll3] - 1;
            double x0 = xstore[inode], y0 = ystore[inode], z0 = zstore[inode];
            if (iz > 0) {
                int idz = (iz - 1) * NGLL * NGLL + iy * NGLL + ix;
                int inodez = ibool[e * NGLL3 + idz] - 1;
                double dz = hypot3(x0 - xstore[inodez], y0 - ystore[inodez], z0 - zstore[inodez]);
                dmin = dmin < dz ? dmin : dz;
            }
            if (iy > 0) {
                int idy = (iz * NGLL + (iy - 1)) * NGLL + ix;
                int inodey = ibool[e * NGLL3 + idy] - 1;
                double dy = hypot3(x0 - xstore[inodey], y0 - ystore[inodey], z0 - zstore[inodey]);
                dmin = dmin < dy ? dmin : dy;
            }
            if (ix > 0) {
                int idx = (iz * NGLL + iy) * NGLL + (ix - 1);
                int inodex = ibool[e * NGLL3 + idx] - 1;
                double dx = hypot3(x0 - xstore[inodex], y0 - ystore[inodex], z0 - zstore[inodex]);
                dmin = dmin < dx ? dmin : dx;
            }
        }}}
    }

    MPI_Allreduce(&dmin, min_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

static uint64_t
get_coordinate_hash(double x, double y, double z,
                    const double *minmax_xyz, double min_dist)
{
    enum { BITS_PER_AXIS = 21 };
    const uint64_t MAX_BINS = UINT64_C(1) << BITS_PER_AXIS;

    const double xmin = minmax_xyz[0], xmax = minmax_xyz[1];
    const double ymin = minmax_xyz[2], ymax = minmax_xyz[3];
    const double zmin = minmax_xyz[4], zmax = minmax_xyz[5];

    const double dx = xmax - xmin;
    const double dy = ymax - ymin;
    const double dz = zmax - zmin;

    if (min_dist <= 0.0 || dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        fprintf(stderr, "Error: invalid input to get_coordinate_hash\n");
        exit(1);
    }

    const double bin_size = min_dist / 10.0;
    uint64_t nbx = (uint64_t)floor(dx / bin_size) + 1;
    uint64_t nby = (uint64_t)floor(dy / bin_size) + 1;
    uint64_t nbz = (uint64_t)floor(dz / bin_size) + 1;

    if (nbx >= MAX_BINS || nby >= MAX_BINS || nbz >= MAX_BINS) {
        fprintf(stderr, "Error: bin count exceeds limit in get_coordinate_hash\n");
        exit(1);
    }

    int64_t xb = llround((x - xmin) / dx * (double)(nbx - 1));
    int64_t yb = llround((y - ymin) / dy * (double)(nby - 1));
    int64_t zb = llround((z - zmin) / dz * (double)(nbz - 1));

    if (xb < 0) xb = 0;
    if (yb < 0) yb = 0;
    if (zb < 0) zb = 0;
    if ((uint64_t)xb >= nbx) xb = (int64_t)(nbx - 1);
    if ((uint64_t)yb >= nby) yb = (int64_t)(nby - 1);
    if ((uint64_t)zb >= nbz) zb = (int64_t)(nbz - 1);

    return ((uint64_t)xb << (2 * BITS_PER_AXIS))
         | ((uint64_t)yb <<      BITS_PER_AXIS)
         |  (uint64_t)zb;
}

/**
 * @brief Build the per-node local-to-global mapping.
 */
static void
create_local2global_mapping(PetscFemCtx *ctx, PetscInt *l2g_out,
                            PetscInt *rstart, PetscInt *rend)
{
    int64_t global_offset = 0;
    int64_t local_owned = (int64_t)ctx->local_owned_dofs;
    MPI_Exscan(&local_owned, &global_offset, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    if (ctx->myrank == 0) global_offset = 0;

    *rstart = (PetscInt)global_offset;
    *rend   = (PetscInt)(global_offset + local_owned);

    int64_t *l2g_indices = (int64_t *)malloc(ctx->nglob * sizeof(int64_t));
    int64_t owned_counter = 0;
    for (int64_t i = 0; i < (int64_t)ctx->nglob; i++) {
        if (ctx->owner_rank[i] == ctx->myrank) {
            l2g_indices[i] = global_offset + owned_counter;
            owned_counter += 1;
        } else {
            l2g_indices[i] = -1;
        }
    }

    if (sizeof(PetscInt) < sizeof(int64_t)) {
        for (int64_t i = 0; i < (int64_t)ctx->nglob; i++) {
            if (l2g_indices[i] > INT_MAX) {
                fprintf(stderr, "Error: global DOF index exceeds PetscInt range; rebuild PETSc with 64-bit indices.\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    int64_t *buf_sd = (int64_t *)malloc(ctx->num_neighbors * ctx->max_nibool * sizeof(int64_t));
    int64_t *buf_rv = (int64_t *)malloc(ctx->num_neighbors * ctx->max_nibool * sizeof(int64_t));
    MPI_Request *req_sd = (MPI_Request *)malloc(ctx->num_neighbors * sizeof(MPI_Request));
    MPI_Request *req_rv = (MPI_Request *)malloc(ctx->num_neighbors * sizeof(MPI_Request));

    for (int i = 0; i < ctx->num_neighbors; i++) {
        int npts = ctx->npts_neighbor[i];
        for (int j = 0; j < npts; j++) {
            int local_node = ctx->send_points[i * ctx->max_nibool + j] - 1;
            buf_sd[i * ctx->max_nibool + j] = l2g_indices[local_node];
        }
        MPI_Isend(buf_sd + i * ctx->max_nibool, npts, MPI_INT64_T,
                  ctx->neighbor_ranks[i], 99, MPI_COMM_WORLD, &req_sd[i]);
        MPI_Irecv(buf_rv + i * ctx->max_nibool, npts, MPI_INT64_T,
                  ctx->neighbor_ranks[i], 99, MPI_COMM_WORLD, &req_rv[i]);
    }

    MPI_Waitall(ctx->num_neighbors, req_rv, MPI_STATUSES_IGNORE);

    for (int i = 0; i < ctx->num_neighbors; i++) {
        for (int j = 0; j < ctx->npts_neighbor[i]; j++) {
            int local_node = ctx->send_points[i * ctx->max_nibool + j] - 1;
            int neighbor_rank = ctx->neighbor_ranks[i];
            if (ctx->owner_rank[local_node] == neighbor_rank) {
                l2g_indices[local_node] = buf_rv[i * ctx->max_nibool + j];
            }
        }
    }
    MPI_Waitall(ctx->num_neighbors, req_sd, MPI_STATUSES_IGNORE);

    for (int i = 0; i < ctx->nglob; i++) {
        l2g_out[i] = (PetscInt)l2g_indices[i];
    }

    free(l2g_indices);
    free(buf_sd); free(buf_rv);
    free(req_sd); free(req_rv);
}

/* =========================================================================
 * COO setup (matrix only)
 * =========================================================================
 * Build the global (i,j) index pattern for every element's dense stiffness
 * block and hand it to PETSc via MatSetPreallocationCOO. The Fortran caller
 * owns the matching values buffer coo_v and passes it to fill_mat_petsc in
 * one shot; that call delegates to MatSetValuesCOO and leaves the matrix
 * in assembled state.
 *
 * The RHS / solution / local-solution vectors and their scatter use the
 * same BlockedLocal pipeline as petsc_routines.c.
 *
 * Index layout for one element follows Fortran's Kloc(dim_out, node_out,
 * dim_in, node_in) column-major linearisation:
 *   k = dim_out + NDIM*node_out + NDIM*NGLL3*dim_in + NDIM*NGLL3*NDIM*node_in
 * with per-element stride dpe*dpe (where dpe = NGLL3*NDIM), so the caller's
 * Kloc(NDIM,NGLL3,NDIM,NGLL3,NSPEC) array can be passed verbatim.
 */
static void
setup_petsc_impl_coo(PetscFemCtx *ctx, const int *ibool)
{
    PetscInt local_owned = (PetscInt)ctx->local_owned_dofs;
    PetscInt rstart, rend;
    PetscInt *l2g_indices = (PetscInt *)malloc(ctx->nglob * sizeof(PetscInt));

    create_local2global_mapping(ctx, l2g_indices, &rstart, &rend);

    int NGLL3 = ctx->NGLL3;
    int NDIM  = ctx->NDIM;
    int dpe   = NGLL3 * NDIM;                            /* dofs per element */
    PetscInt local_dofs = local_owned * NDIM;

    PetscCount coo_n = (PetscCount)ctx->nspec * (PetscCount)dpe * (PetscCount)dpe;
    PetscInt *coo_i  = (PetscInt *)malloc((size_t)(coo_n > 0 ? coo_n : 1) * sizeof(PetscInt));
    PetscInt *coo_j  = (PetscInt *)malloc((size_t)(coo_n > 0 ? coo_n : 1) * sizeof(PetscInt));

    /* Matrix COO indices. PETSc accepts off-rank rows/cols in the COO arrays
     * and routes the corresponding values to their owners at assembly time. */
    for (int e = 0; e < ctx->nspec; e++) {
        PetscCount e_off = (PetscCount)e * (PetscCount)dpe * (PetscCount)dpe;
        for (int node_in = 0; node_in < NGLL3; node_in++) {
            PetscInt gid_in = l2g_indices[ibool[e * NGLL3 + node_in] - 1];
            for (int dim_in = 0; dim_in < NDIM; dim_in++) {
                PetscInt jcol = gid_in * NDIM + dim_in;
                for (int node_out = 0; node_out < NGLL3; node_out++) {
                    PetscInt gid_out = l2g_indices[ibool[e * NGLL3 + node_out] - 1];
                    for (int dim_out = 0; dim_out < NDIM; dim_out++) {
                        PetscInt irow = gid_out * NDIM + dim_out;
                        PetscCount k_local = (PetscCount)dim_out
                                           + (PetscCount)NDIM * node_out
                                           + (PetscCount)NDIM * NGLL3 * dim_in
                                           + (PetscCount)NDIM * NGLL3 * NDIM * node_in;
                        PetscCount k = e_off + k_local;
                        coo_i[k] = irow;
                        coo_j[k] = jcol;
                    }
                }
            }
        }
    }

    /* Create matrix; default to MATAIJ since COO maps naturally onto scalar
     * AIJ. -mat_type aijcusparse / aijkokkos / ... still works through
     * MatSetFromOptions. */
    PetscCallVoid(MatCreate(PETSC_COMM_WORLD, &ctx->K_mat));
    PetscCallVoid(MatSetSizes(ctx->K_mat, local_dofs, local_dofs, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCallVoid(MatSetBlockSize(ctx->K_mat, NDIM));
    PetscCallVoid(MatSetType(ctx->K_mat, MATAIJ));
    PetscCallVoid(MatSetFromOptions(ctx->K_mat));
    PetscCallVoid(MatSetPreallocationCOO(ctx->K_mat, coo_n, coo_i, coo_j));

    /* coo_i / coo_j have been absorbed into PETSc's internal permutation
     * table; the values array is supplied later by the Fortran caller. */
    free(coo_i); free(coo_j);

    /* Node-level local-to-global map for the vector BlockedLocal API and
     * the global-to-local scatter. */
    ISLocalToGlobalMapping l2g_mapping;
    PetscCallVoid(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, NDIM, ctx->nglob,
                                               l2g_indices, PETSC_COPY_VALUES, &l2g_mapping));

    /* RHS / solution vectors — original BlockedLocal pipeline. */
    PetscCallVoid(VecCreate(PETSC_COMM_WORLD, &ctx->rhs_vec));
    PetscCallVoid(VecSetSizes(ctx->rhs_vec, local_dofs, PETSC_DETERMINE));
    PetscCallVoid(VecSetBlockSize(ctx->rhs_vec, NDIM));
    PetscCallVoid(VecSetType(ctx->rhs_vec, VECMPI));
    PetscCallVoid(VecSetFromOptions(ctx->rhs_vec));
    PetscCallVoid(VecSetUp(ctx->rhs_vec));
    PetscCallVoid(VecSetLocalToGlobalMapping(ctx->rhs_vec, l2g_mapping));

    PetscCallVoid(VecDuplicate(ctx->rhs_vec, &ctx->sol_vec));

    /* Local solution vector and scatter (unchanged from AIJ version). */
    PetscCallVoid(VecCreate(PETSC_COMM_SELF, &ctx->local_sol_vec));
    PetscCallVoid(VecSetSizes(ctx->local_sol_vec, ctx->nglob * NDIM, PETSC_DETERMINE));
    PetscCallVoid(VecSetBlockSize(ctx->local_sol_vec, NDIM));
    PetscCallVoid(VecSetType(ctx->local_sol_vec, VECSEQ));
    PetscCallVoid(VecSetFromOptions(ctx->local_sol_vec));
    PetscCallVoid(VecSetUp(ctx->local_sol_vec));

    IS from_is, to_is;
    PetscCallVoid(ISCreateBlock(PETSC_COMM_SELF, NDIM, ctx->nglob,
                                l2g_indices, PETSC_COPY_VALUES, &from_is));
    PetscInt *seq_indices = (PetscInt *)malloc(ctx->nglob * sizeof(PetscInt));
    for (int i = 0; i < ctx->nglob; i++) seq_indices[i] = i;
    PetscCallVoid(ISCreateBlock(PETSC_COMM_SELF, NDIM, ctx->nglob,
                                seq_indices, PETSC_COPY_VALUES, &to_is));
    PetscCallVoid(VecScatterCreate(ctx->sol_vec, from_is, ctx->local_sol_vec, to_is, &ctx->sol_scatter));

    PetscCallVoid(ISDestroy(&from_is));
    PetscCallVoid(ISDestroy(&to_is));
    PetscCallVoid(ISLocalToGlobalMappingDestroy(&l2g_mapping));
    free(seq_indices);
    free(l2g_indices);

    // allocate connectivity matrix
    ctx->conn = (PetscInt *)malloc((size_t)(ctx->nspec * ctx->NGLL3 * sizeof(PetscInt)));
    for(size_t i = 0; i < (size_t)ctx->nspec * ctx->NGLL3; i++) {
        ctx->conn[i] = ibool[i] - 1; // convert to 0-based indexing for C
    }
}

/**
 * @brief Fortran entry point: build ownership info and create PETSc objects.
 */
void FC_FUNC_(setup_petsc, SETUP_PETSC)(long *h,
                                        int nglob, int myrank,
                                        int nspec, int NGLL3, int NDIM,
                                        int num_interfaces_ext_mesh,
                                        const int *nibool_interfaces_ext_mesh,
                                        const int *ibool_interfaces_ext_mesh,
                                        int max_nibool,
                                        const int *my_neighbors_ext_mesh,
                                        const float *xstore,
                                        const float *ystore,
                                        const float *zstore,
                                        const int *ibool,
                                        int *owner_rank)
{
    PetscFemCtx *ctx = (PetscFemCtx *)calloc(1, sizeof(PetscFemCtx));
    *h = (long)ctx;

    if (myrank == 0) {
        printf("-------------------------------\n");
        printf("creating ownerships (COO) ...\n");
        printf("-------------------------------\n");
    }

    for (int i = 0; i < nglob; i++) owner_rank[i] = myrank;

    double minmax_xyz[6], min_dist;
    get_xyz_minmax(nspec, NGLL3, nglob, ibool, xstore, ystore, zstore, minmax_xyz, &min_dist);

    ctx->num_neighbors  = num_interfaces_ext_mesh;
    ctx->max_nibool     = max_nibool;
    ctx->neighbor_ranks = my_neighbors_ext_mesh;
    ctx->npts_neighbor  = nibool_interfaces_ext_mesh;
    ctx->send_points    = ibool_interfaces_ext_mesh;

    /* Lowest-score-wins tie-break on interface nodes. */
    for (int ineigh = 0; ineigh < ctx->num_neighbors; ineigh++) {
        int neigh_rank = ctx->neighbor_ranks[ineigh];
        int num_points = ctx->npts_neighbor[ineigh];
        const int *points = &ctx->send_points[ineigh * ctx->max_nibool];

        for (int ipoin = 0; ipoin < num_points; ipoin++) {
            int node = points[ipoin] - 1;
            double x = (double)xstore[node];
            double y = (double)ystore[node];
            double z = (double)zstore[node];
            uint64_t hash_val   = get_coordinate_hash(x, y, z, minmax_xyz, min_dist);
            uint64_t my_score   = get_hash_score(hash_val, (uint64_t)owner_rank[node]);
            uint64_t neigh_score = get_hash_score(hash_val, (uint64_t)neigh_rank);
            if (neigh_score > my_score) {
                owner_rank[node] = neigh_rank;
            }
        }
    }

    ctx->local_owned_dofs = 0;
    for (int i = 0; i < nglob; i++) {
        if (owner_rank[i] == myrank) ctx->local_owned_dofs++;
    }

    printf("Rank %d (COO) owns %d / %d nodes\n",
           myrank, (int)ctx->local_owned_dofs, nglob);

    ctx->nglob      = nglob;
    ctx->owner_rank = owner_rank;
    ctx->myrank     = myrank;
    ctx->nspec      = nspec;
    ctx->NGLL3      = NGLL3;
    ctx->NDIM       = NDIM;

    setup_petsc_impl_coo(ctx, ibool);
}

/**
 * @brief Push the caller's full COO values buffer into the global stiffness
 *        matrix.
 *
 * @p coo_v must have length nspec * (NGLL3*NDIM)^2 and follow the same
 * per-element layout that setup_petsc_impl_coo used when building the COO
 * (i,j) pattern — namely, Fortran column-major Kloc(NDIM,NGLL3,NDIM,NGLL3)
 * per element, contiguous across elements.
 *
 * MatSetValuesCOO sums duplicate (i,j) entries, routes off-rank rows/cols
 * to their owners, and leaves the matrix in assembled state, so no further
 * MatAssemblyBegin/End is needed.
 */
void FC_FUNC_(fill_mat_petsc, FILL_MAT_PETSC)(long *h,
                                              const PetscScalar *coo_v)
{
    PetscFemCtx *ctx = CTX(h);
    PetscCallVoid(MatSetValuesCOO(ctx->K_mat, coo_v, ADD_VALUES));
}

/**
 * @brief Accumulate element-wise RHS contributions into the global RHS vector. The caller's f_elem
 *        VecSetValuesBlockedLocal + ADD_VALUES (unchanged from petsc_routines.c).
 * @param f_elem with shape(NSPEC, NGLL3, NDIM) per element
 */
void FC_FUNC_(fill_vec_petsc, FILL_VEC_PETSC)(long *h,
                                              const PetscScalar *f_elem)
{
    PetscFemCtx *ctx = CTX(h);
    int NGLL3 = ctx->NGLL3;
    int nspec = ctx->nspec;

    PetscCallVoid(VecSetValuesBlockedLocal(ctx->rhs_vec,
                                           NGLL3*nspec, ctx->conn,
                                           f_elem,
                                           ADD_VALUES));

    PetscCallVoid(VecAssemblyBegin(ctx->rhs_vec));
    PetscCallVoid(VecAssemblyEnd  (ctx->rhs_vec));
}

/**
 * @brief Solve K u = f via PETSc KSP.
 */
void FC_FUNC_(solve_petsc, SOLVE_PETSC)(long *h)
{
    PetscFemCtx *ctx = CTX(h);

    PetscCallVoid(KSPCreate(PETSC_COMM_WORLD, &ctx->ksp_solver));
    PetscCallVoid(KSPSetOperators(ctx->ksp_solver, ctx->K_mat, ctx->K_mat));
    PetscCallVoid(KSPSetFromOptions(ctx->ksp_solver));
    PetscCallVoid(KSPSolve(ctx->ksp_solver, ctx->rhs_vec, ctx->sol_vec));
}

/**
 * @brief Copy the global solution into the Fortran @p displ array, populating
 *        ghost nodes via the PETSc scatter built in setup_petsc_impl_coo.
 */
void FC_FUNC_(extract_petsc, EXTRACT_PETSC)(long *h, double *displ)
{
    PetscFemCtx *ctx = CTX(h);
    int NDIM = ctx->NDIM;

    PetscCallVoid(VecScatterBegin(ctx->sol_scatter, ctx->sol_vec, ctx->local_sol_vec, INSERT_VALUES, SCATTER_FORWARD));
    PetscCallVoid(VecScatterEnd  (ctx->sol_scatter, ctx->sol_vec, ctx->local_sol_vec, INSERT_VALUES, SCATTER_FORWARD));

    const PetscScalar *array;
    PetscCallVoid(VecGetArrayRead(ctx->local_sol_vec, &array));
    for (int i = 0; i < ctx->nglob * NDIM; i++) {
        displ[i] = (double)array[i];
    }
    PetscCallVoid(VecRestoreArrayRead(ctx->local_sol_vec, &array));
}

/**
 * @brief Destroy PETSc objects, free COO buffers, null out the handle.
 */
void FC_FUNC_(cleanup_petsc, CLEANUP_PETSC)(long *h)
{
    PetscFemCtx *ctx = CTX(h);

    PetscCallVoid(KSPDestroy(&ctx->ksp_solver));
    PetscCallVoid(VecDestroy(&ctx->rhs_vec));
    PetscCallVoid(VecDestroy(&ctx->sol_vec));
    PetscCallVoid(MatDestroy(&ctx->K_mat));
    PetscCallVoid(VecDestroy(&ctx->local_sol_vec));
    PetscCallVoid(VecScatterDestroy(&ctx->sol_scatter));

    free(ctx->conn);
    free(ctx);
    *h = 0;
}
