# CMakeLists.txt for static solver
#
# Builds xstatic3D executable.
# Reuses all specfem3D solver objects (via specfem3d_solver_objects)
# except the main program specfem3D.f90.

###############################################################################
# Executable: xstatic3D
###############################################################################

add_executable(xstatic3D
  petsc_interfaces.f90
  static_module.f90
  xstatic3D.F90
)

if(WITH_PETSC)
  target_sources(xstatic3D PRIVATE petsc_routines.c)
  target_link_libraries(xstatic3D PRIVATE PETSc::PETSc)
else()
  target_sources(xstatic3D PRIVATE petsc_routines_stubs.c)
endif()

target_link_libraries(xstatic3D PRIVATE
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(WITH_MPI)
  target_link_libraries(xstatic3D PRIVATE MPI::MPI_Fortran)
endif()

if(HAS_GPU)
  target_link_libraries(xstatic3D PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_VTK)
  target_link_libraries(xstatic3D PRIVATE ${VTK_LIBRARIES})
endif()

# Install
install(TARGETS xstatic3D DESTINATION bin)
