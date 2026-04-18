program prog_static_centrifugal
  use static_module
  use specfem_par

#ifdef WITH_PETSC
  use petscsys
#endif

  implicit none

  integer :: ierr

  ! MPI initialization
  call init_mpi()

#ifdef WITH_PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  ! reads in parameters
  call initialize_simulation()

  ! reads in external mesh
  call read_mesh_databases()

  ! reads in moho mesh
  call read_mesh_databases_moho()

  ! reads adjoint parameters
  call read_mesh_databases_adjoint()

  ! sets up local time stepping (LTS)
  if (LTS_MODE) call lts_setup()

  ! for coupling with external codes
  call couple_with_injection_setup()

  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! detects surfaces
  call detect_mesh_surfaces()

  ! prepares sources and receivers
  call setup_sources_receivers()

  ! sets up and precomputes simulation arrays
  call prepare_timerun()

  ! run static solver
  call initialize_static_module()
  call solve_static_problem()
  call finalize_static_module()

  ! saves last time frame and finishes kernel calculations
  call finalize_simulation()

#ifdef WITH_PETSC
  call PetscFinalize(ierr)
#endif

  ! MPI finish
  call finalize_mpi()


  
end program prog_static_centrifugal