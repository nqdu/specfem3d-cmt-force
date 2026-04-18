# CMakeLists.txt for specfem3d_shared library
#
# Builds shared utility code used by multiple SPECFEM3D executables.

###############################################################################
# Always-included Fortran sources
###############################################################################

# shared_par.f90 is a MODULE that all other sources depend on – listed first
# so CMake's Fortran dependency scanner picks up the ordering.
set(SHARED_FORTRAN_SOURCES
  shared_par.F90
  assemble_MPI_scalar.f90
  check_mesh_resolution.f90
  count_number_of_sources.f90
  create_name_database.f90
  define_derivation_matrices.f90
  define_mass_matrices.f90
  detect_surface.f90
  exit_mpi.f90
  fft.f90
  get_attenuation_model.f90
  get_element_face.f90
  get_global.f90
  get_jacobian_boundaries.f90
  get_shape2D.f90
  get_shape3D.f90
  gll_library.f90
  hex_nodes.f90
  heap_sort.f90
  init_openmp.F90
  lagrange_poly.f90
  merge_sort.f90
  netlib_specfun_erf.f90
  prepare_assemble_MPI.f90
  read_parameter_file.F90
  read_topo_bathy_file.f90
  read_value_parameters.f90
  recompute_jacobian.f90
  safe_alloc_mod.f90
  save_header_file.f90
  search_kdtree.f90
  sort_array_coordinates.f90
  utm_geo.f90
  write_VTK_data.f90
)

# Module files that depend on shared_par (compiled with non-MPI compiler)
list(APPEND SHARED_FORTRAN_SOURCES
  adios_manager.F90
  hdf5_manager.F90
)

###############################################################################
# Always-included C sources
###############################################################################

set(SHARED_C_SOURCES
  force_ftz.c
  param_reader.c
  write_c_binary.c
)

###############################################################################
# Conditional: MPI vs serial
###############################################################################

if(WITH_MPI)
  list(APPEND SHARED_FORTRAN_SOURCES parallel.f90)
else()
  list(APPEND SHARED_FORTRAN_SOURCES serial.f90)
endif()

###############################################################################
# Conditional: ADIOS / ADIOS2
###############################################################################

if(WITH_ADIOS OR WITH_ADIOS2)
  list(APPEND SHARED_C_SOURCES adios_helpers_addons.c)
  list(APPEND SHARED_FORTRAN_SOURCES
    adios_helpers_definitions.F90
    adios_helpers_readers.F90
    adios_helpers_writers.F90
    adios_helpers.f90
  )
else()
  list(APPEND SHARED_C_SOURCES adios_method_stubs.c)
endif()

###############################################################################
# Conditional: ASDF
###############################################################################

if(WITH_ASDF)
  list(APPEND SHARED_FORTRAN_SOURCES asdf_manager.f90)
else()
  list(APPEND SHARED_C_SOURCES asdf_method_stubs.c)
  list(APPEND SHARED_FORTRAN_SOURCES asdf_manager_stubs.f90)
endif()

###############################################################################
# Library target
###############################################################################

add_library(specfem3d_shared STATIC
  ${SHARED_FORTRAN_SOURCES}
  ${SHARED_C_SOURCES}
)

###############################################################################
# Include directories
###############################################################################

target_include_directories(specfem3d_shared PUBLIC
  ${CMAKE_BINARY_DIR}/setup
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

###############################################################################
# MPI
###############################################################################

if(WITH_MPI)
  target_compile_definitions(specfem3d_shared PUBLIC WITH_MPI)
  target_link_libraries(specfem3d_shared PUBLIC MPI::MPI_Fortran MPI::MPI_C)
endif()

###############################################################################
# Vectorization
###############################################################################

if(FORCE_VECTORIZATION)
  target_compile_definitions(specfem3d_shared PUBLIC FORCE_VECTORIZATION)
endif()

###############################################################################
# OpenMP
###############################################################################

if(WITH_OPENMP)
  target_link_libraries(specfem3d_shared PUBLIC OpenMP::OpenMP_Fortran)
  target_compile_definitions(specfem3d_shared PUBLIC USE_OPENMP)
endif()

###############################################################################
# HDF5
###############################################################################

if(WITH_HDF5)
  target_compile_definitions(specfem3d_shared PUBLIC USE_HDF5)
  target_link_libraries(specfem3d_shared PUBLIC ${HDF5_Fortran_LIBRARIES})
  target_include_directories(specfem3d_shared PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()

###############################################################################
# ADIOS / ADIOS2
###############################################################################

if(WITH_ADIOS OR WITH_ADIOS2)
  target_include_directories(specfem3d_shared PUBLIC
    ${ADIOS_INCLUDE_DIRS}
    ${ADIOS2_INCLUDE_DIRS}
  )
  if(WITH_ADIOS)
    target_compile_definitions(specfem3d_shared PUBLIC USE_ADIOS)
  endif()
  if(WITH_ADIOS2)
    target_compile_definitions(specfem3d_shared PUBLIC USE_ADIOS2)
  endif()
endif()
