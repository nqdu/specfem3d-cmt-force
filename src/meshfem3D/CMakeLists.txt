###############################################################################
# meshfem3D - builds xmeshfem3D executable
###############################################################################

# Module files (must be compiled first for module dependencies)
set(MESHFEM3D_MODULE_SOURCES
  meshfem3D_par.f90
  chunk_earth_mesh_mod.f90
)

# Main source files
set(MESHFEM3D_SOURCES
  meshfem3D.F90
  calc_gll_points.f90
  check_mesh_quality.f90
  compute_parameters.f90
  create_meshfem_mesh.f90
  create_CPML_regions.f90
  create_interfaces_mesh.f90
  create_visual_files.f90
  define_subregions.f90
  define_subregions_heuristic.f90
  define_superbrick.f90
  determine_cavity.f90
  earth_chunk.f90
  get_flags_boundaries.f90
  get_MPI_cutplanes_eta.f90
  get_MPI_cutplanes_xi.f90
  get_wavefield_discontinuity.f90
  read_mesh_parameter_file.f90
  read_value_mesh_parameters.f90
  save_databases.F90
  save_mesh_as_cubit.f90
  store_boundaries.f90
  store_coords.f90
)

# HDF5 source (always compiled, uses preprocessor for conditional HDF5 support)
set(MESHFEM3D_HDF5_SOURCES
  save_databases_hdf5.F90
)

# ADIOS sources (conditional)
set(MESHFEM3D_ADIOS_SOURCES
  save_databases_adios.F90
)

# Combine all sources
set(MESHFEM3D_ALL_SOURCES
  ${MESHFEM3D_MODULE_SOURCES}
  ${MESHFEM3D_SOURCES}
  ${MESHFEM3D_HDF5_SOURCES}
)

if(WITH_ADIOS OR WITH_ADIOS2)
  list(APPEND MESHFEM3D_ALL_SOURCES ${MESHFEM3D_ADIOS_SOURCES})
endif()

add_executable(xmeshfem3D ${MESHFEM3D_ALL_SOURCES})

target_link_libraries(xmeshfem3D PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xmeshfem3D PRIVATE MPI::MPI_Fortran)
endif()

if(WITH_HDF5)
  target_link_libraries(xmeshfem3D PRIVATE ${HDF5_Fortran_LIBRARIES})
  target_include_directories(xmeshfem3D PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
endif()

if(WITH_ADIOS)
  target_link_libraries(xmeshfem3D PRIVATE ${ADIOS_LIBRARIES})
  target_include_directories(xmeshfem3D PRIVATE ${ADIOS_INCLUDE_DIRS})
elseif(WITH_ADIOS2)
  target_link_libraries(xmeshfem3D PRIVATE ${ADIOS2_LIBRARIES})
  target_include_directories(xmeshfem3D PRIVATE ${ADIOS2_INCLUDE_DIRS})
endif()

# Install
install(TARGETS xmeshfem3D DESTINATION bin)
