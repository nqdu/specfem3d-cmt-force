!> @file petsc_interfaces.f90
!> @brief Fortran iso_c_binding interfaces to the C PETSc FEM routines.
!>
!> Two interchangeable C implementations expose the symbols declared here:
!>   - petsc_routines.c     : AIJ/BAIJ MatSetValuesBlockedLocal assembly
!>   - petsc_routines_coo.c : MatSetPreallocationCOO / MatSetValuesCOO
!>                            assembly for the stiffness matrix (the RHS
!>                            vector still uses VecSetValuesBlockedLocal)
!>
!> Both implementations share the same Fortran-callable API, so callers
!> never need to know which one is linked.
!>
!> Usage pattern:
!> @code
!>   use petsc_interfaces
!>   integer(c_long) :: h
!>   call setup_petsc(...)
!>   ...
!>   call cleanup_petsc(h)
!> @endcode

module petsc_interfaces
  use iso_c_binding, only: c_long, c_int, c_float, c_double
  implicit none

  interface

    !> @brief Build ownership metadata and create/preallocate PETSc objects.
    !>
    !> @param[out] h                           Solver context handle.
    !> @param[in]  nglob                        Total nodes on this partition.
    !> @param[in]  myrank                       MPI rank of this process.
    !> @param[in]  nspec                        Number of local spectral elements.
    !> @param[in]  NGLL3                        Number of GLL points per element.
    !> @param[in]  NDIM                         Number of spatial dimensions.
    !> @param[in]  num_interfaces_ext_mesh      Number of shared interfaces.
    !> @param[in]  nibool_interfaces_ext_mesh   Node count per interface.
    !> @param[in]  ibool_interfaces_ext_mesh    1-based node indices (column-major).
    !> @param[in]  max_nibool                   Leading dimension of ibool array.
    !> @param[in]  my_neighbors_ext_mesh        MPI rank of each interface neighbour.
    !> @param[in]  xstore                       X coordinates of local/global nodes.
    !> @param[in]  ystore                       Y coordinates of local/global nodes.
    !> @param[in]  zstore                       Z coordinates of local/global nodes.
    !> @param[in]  ibool                        Element connectivity (1-based, column-major).
    !> @param[out] owner_rank                   Owning rank per node (length nglob).
    subroutine setup_petsc(h, nglob, myrank, &
                           nspec, NGLL3, NDIM, &
                           num_interfaces_ext_mesh, &
                           nibool_interfaces_ext_mesh, &
                           ibool_interfaces_ext_mesh, &
                           max_nibool, &
                           my_neighbors_ext_mesh, &
                           xstore, ystore, zstore, &
                           ibool, owner_rank) bind(C, name="setup_petsc_")
      import c_long, c_int, c_float
      implicit none
      integer(c_long), intent(out)       :: h
      integer(c_int),  intent(in), value :: nglob
      integer(c_int),  intent(in), value :: myrank
      integer(c_int),  intent(in), value :: nspec
      integer(c_int),  intent(in), value :: NGLL3
      integer(c_int),  intent(in), value :: NDIM
      integer(c_int),  intent(in), value :: num_interfaces_ext_mesh
      integer(c_int),  intent(in)        :: nibool_interfaces_ext_mesh(*)
      integer(c_int),  intent(in)        :: ibool_interfaces_ext_mesh(*)
      integer(c_int),  intent(in), value :: max_nibool
      integer(c_int),  intent(in)        :: my_neighbors_ext_mesh(*)
      real(c_float),   intent(in)        :: xstore(*)
      real(c_float),   intent(in)        :: ystore(*)
      real(c_float),   intent(in)        :: zstore(*)
      integer(c_int),  intent(in)        :: ibool(*)
      integer(c_int),  intent(out)       :: owner_rank(*)
    end subroutine setup_petsc

    !> @brief Push the full element-wise stiffness contribution into K.
    !>
    !> COO backend: @p coo_v has length NSPEC*(NGLL3*NDIM)**2 and is laid out
    !> per element as Fortran Kloc(NDIM,NGLL3,NDIM,NGLL3,NSPEC) (i.e. exactly
    !> the layout produced by writing each element's Kloc block contiguously
    !> in ispec order). The call delegates to MatSetValuesCOO, which sums
    !> duplicate (i,j) entries, routes off-rank values, and leaves K
    !> assembled.
    !>
    !> AIJ backend (petsc_routines.c): the same @p coo_v buffer is used in a
    !> per-block loop through MatSetValuesBlockedLocal — see that file for
    !> details. Stubs (no-PETSc build) ignore the argument.
    !>
    !> @param[in] h     Solver context handle.
    !> @param[in] coo_v Stiffness values, length NSPEC*(NGLL3*NDIM)**2.
    subroutine fill_mat_petsc(h, coo_v) &
        bind(C, name="fill_mat_petsc_")
      import c_long, c_double
      implicit none
      integer(c_long), intent(in) :: h
      real(c_double),  intent(in) :: coo_v(*)
    end subroutine fill_mat_petsc

    !> @brief Accumulate one element's force vector into the global RHS vector.
    !>
    !> Identical in both backends: values are added via VecSetValuesBlockedLocal
    !> with ADD_VALUES and merged across ranks inside assemble_petsc.
    !>
    !> @param[in] h             Solver context handle.
    !> @param[in] f_elem        Element force values (length NSPEC*NGLL3*NDIM, Fortran order: f_elem(NDIM,NGLL3,NSPEC) i.e. ispec-fastest, then GLL, then dim
    !>                          DOF-fastest order).
    subroutine fill_vec_petsc(h, f_elem) &
        bind(C, name="fill_vec_petsc_")
      import c_long, c_int, c_double
      implicit none
      integer(c_long), intent(in) :: h
      real(c_double),  intent(in) :: f_elem(*)
    end subroutine fill_vec_petsc

    !> @brief Solve the assembled linear system K*u = f.
    !> @param[in] h  Solver context handle.
    subroutine solve_petsc(h) bind(C, name="solve_petsc_")
      import c_long
      implicit none
      integer(c_long), intent(in) :: h
    end subroutine solve_petsc

    !> @brief Extract solution from PETSc and scatter to ghost nodes via MPI.
    !>
    !> @param[in]     h      Solver context handle.
    !> @param[in,out] displ  Native displacement array (size nglob*NDIM); owned entries
    !>                       are filled from PETSc then ghost entries are updated via MPI.
    subroutine extract_petsc(h, displ) bind(C, name="extract_petsc_")
      import c_long, c_double
      implicit none
      integer(c_long), intent(in)    :: h
      real(c_double),  intent(inout) :: displ(*)
    end subroutine extract_petsc

    !> @brief Destroy all PETSc objects, finalise the library, and free the context.
    !> @param[in,out] h  Solver context handle; set to 0 on return.
    subroutine cleanup_petsc(h) bind(C, name="cleanup_petsc_")
      import c_long
      implicit none
      integer(c_long), intent(inout) :: h
    end subroutine cleanup_petsc

  end interface

end module petsc_interfaces
