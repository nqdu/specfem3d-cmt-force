!> @file petsc_interfaces.f90
!> @brief Fortran iso_c_binding interfaces to the C PETSc FEM routines in petsc_routines.c.
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

    !> @brief Accumulate one element's stiffness block into the global matrix.
    !>
    !> @param[in] h              Solver context handle.
    !> @param[in] global_indices 0-based global node indices (length NGLL3).
    !> @param[in] k_elem         Dense element stiffness values (row-major, NGLL3*NDIM x NGLL3*NDIM).
    subroutine fill_mat_petsc(h, global_indices, k_elem) &
        bind(C, name="fill_mat_petsc_")
      import c_long, c_int, c_double
      implicit none
      integer(c_long), intent(in) :: h
      integer(c_int),  intent(in) :: global_indices(*)
      real(c_double),  intent(in) :: k_elem(*)
    end subroutine fill_mat_petsc

    !> @brief Accumulate one element's force vector into the global RHS vector.
    !>
    !> @param[in] h              Solver context handle.
    !> @param[in] global_indices 0-based global node indices (length NGLL3).
    !> @param[in] f_elem         Element force values (length NGLL3*NDIM, DOF-then-node order).
    subroutine fill_vec_petsc(h, global_indices, f_elem) &
        bind(C, name="fill_vec_petsc_")
      import c_long, c_int, c_double
      implicit none
      integer(c_long), intent(in) :: h
      integer(c_int),  intent(in) :: global_indices(*)
      real(c_double),  intent(in) :: f_elem(*)
    end subroutine fill_vec_petsc

    !> @brief Finalise global matrix and RHS vector assembly (triggers MPI communication).
    !> @param[in] h  Solver context handle.
    subroutine assemble_petsc(h) bind(C, name="assemble_petsc_")
      import c_long
      implicit none
      integer(c_long), intent(in) :: h
    end subroutine assemble_petsc

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
