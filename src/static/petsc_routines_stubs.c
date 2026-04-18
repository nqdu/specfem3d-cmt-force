/**
 * @file petsc_routines_nopetsc.c
 * @brief No-op stubs for all PETSc interface routines.
 *
 * Compiled instead of petsc_routines.c when PETSc is not available
 * (i.e. when WITH_PETSC is not set).  Every function is a silent no-op
 * so the rest of the code links and runs without modification; callers
 * that actually invoke these routines will simply do nothing.
 */

#include <stddef.h>   /* NULL */
#include <stdio.h>    /* fprintf, stderr */
#include "config.h"

/* Suppress unused-parameter warnings for all stub arguments. */
#define UNUSED(x) (void)(x)

/**
 * @brief No-op stub: does not create any matrix, vectors, or ownership metadata.
 */
void FC_FUNC_(setup_petsc,SETUP_PETSC)(long *h,
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
    UNUSED(h); UNUSED(nglob); UNUSED(myrank);
    UNUSED(nspec); UNUSED(NGLL3); UNUSED(NDIM);
    UNUSED(num_interfaces_ext_mesh); UNUSED(nibool_interfaces_ext_mesh);
    UNUSED(ibool_interfaces_ext_mesh); UNUSED(max_nibool);
    UNUSED(my_neighbors_ext_mesh);
    UNUSED(xstore); UNUSED(ystore); UNUSED(zstore);
    UNUSED(ibool); UNUSED(owner_rank);
}

/**
 * @brief No-op stub: does not assemble any element stiffness contributions.
 */
void FC_FUNC_(fill_mat_petsc,FILL_MAT_PETSC)(long *h,
                                             const int    *global_indices,
                                             const double *k_elem)
{
    UNUSED(h); UNUSED(global_indices); UNUSED(k_elem);
}

/**
 * @brief No-op stub: does not assemble any element force contributions.
 */
void FC_FUNC_(fill_vec_petsc,FILL_VEC_PETSC)(long *h,
                                             const int    *global_indices,
                                             const double *f_elem)
{
    UNUSED(h); UNUSED(global_indices); UNUSED(f_elem);
}

/**
 * @brief No-op stub: does not finalise any assembly.
 */
void FC_FUNC_(assemble_petsc,ASSEMBLE_PETSC)(long *h)
{
    UNUSED(h);
}

/**
 * @brief No-op stub: does not solve any system.
 */
void FC_FUNC_(solve_petsc,SOLVE_PETSC)(long *h)
{
    UNUSED(h);
}

/**
 * @brief No-op stub: does not extract any data.
 */
void FC_FUNC_(extract_petsc,EXTRACT_PETSC)(long *h, double *displ)
{
    UNUSED(h); UNUSED(displ);
}

/**
 * @brief No-op stub: sets handle to 0.
 * @param[in,out] h  Set to 0 on return.
 */
void FC_FUNC_(cleanup_petsc,CLEANUP_PETSC)(long *h)
{
    *h = 0;
}
