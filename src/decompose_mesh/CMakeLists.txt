###############################################################################
# decompose_mesh - builds xdecompose_mesh and xdecompose_mesh_mpi executables
###############################################################################

# Module files (must be compiled first for module dependencies)
set(DECOMPOSE_MODULE_SOURCES
  decompose_mesh_par.F90
  fault_scotch.f90
  part_decompose_mesh.F90
)

# Main source files for xdecompose_mesh (serial version)
set(DECOMPOSE_MESH_SOURCES
  check_valence.F90
  decompose_mesh.F90
  lts_helper.F90
  lts_setup_elements.F90
  partition_scotch.F90
  partition_metis.F90
  partition_patoh.F90
  partition_rows.F90
  program_decompose_mesh.f90
  read_mesh_files.F90
  write_mesh_databases.F90
)

# C/C++ sources
set(DECOMPOSE_C_SOURCES
  wrap_metis.c
)
set(DECOMPOSE_CXX_SOURCES
  wrap_patoh.cpp
)

# HDF5 sources (always compiled, uses preprocessor for conditional HDF5 support)
set(DECOMPOSE_HDF5_SOURCES
  part_decompose_mesh_hdf5.F90
  write_mesh_databases_hdf5.F90
)

# MPI version sources
set(DECOMPOSE_MPI_MODULE_SOURCES
  module_qsort.f90
  module_database.f90
  module_mesh.f90
  module_partition.f90
)

set(DECOMPOSE_MPI_SOURCES
  program_decompose_mesh_mpi.f90
)

###############################################################################
# Shared module object library (compiled once to prevent parallel .mod races)
###############################################################################

add_library(decompose_mesh_par_obj OBJECT
  ${DECOMPOSE_MODULE_SOURCES}
  ${DECOMPOSE_HDF5_SOURCES}
)
# Must wait for specfem3d_shared so that constants.mod and other modules exist
add_dependencies(decompose_mesh_par_obj specfem3d_shared)
target_link_libraries(decompose_mesh_par_obj PRIVATE specfem3d_shared)
target_include_directories(decompose_mesh_par_obj PRIVATE
  ${CMAKE_BINARY_DIR}/setup
  ${CMAKE_BINARY_DIR}/obj)
if(WITH_SCOTCH)
  target_include_directories(decompose_mesh_par_obj PRIVATE ${SCOTCH_INCLUDE_DIRS})
  target_compile_definitions(decompose_mesh_par_obj PRIVATE USE_SCOTCH)
endif()
if(WITH_MPI)
  target_link_libraries(decompose_mesh_par_obj PRIVATE MPI::MPI_Fortran)
endif()
if(WITH_HDF5)
  target_include_directories(decompose_mesh_par_obj PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
endif()

###############################################################################
# xdecompose_mesh (serial version - always built)
###############################################################################

add_executable(xdecompose_mesh
  $<TARGET_OBJECTS:decompose_mesh_par_obj>
  ${DECOMPOSE_MESH_SOURCES}
  ${DECOMPOSE_C_SOURCES}
  ${DECOMPOSE_CXX_SOURCES}
)

target_link_libraries(xdecompose_mesh PRIVATE specfem3d_shared)

# Scotch
if(WITH_SCOTCH)
  target_include_directories(xdecompose_mesh PRIVATE ${SCOTCH_INCLUDE_DIRS})
  target_link_libraries(xdecompose_mesh PRIVATE ${SCOTCH_LIBRARIES})
  target_compile_definitions(xdecompose_mesh PRIVATE USE_SCOTCH)
endif()

# METIS
if(WITH_METIS)
  target_include_directories(xdecompose_mesh PRIVATE ${METIS_INCLUDE_DIRS})
  target_link_libraries(xdecompose_mesh PRIVATE ${METIS_LIBRARIES})
  target_compile_definitions(xdecompose_mesh PRIVATE USE_METIS)
endif()

# PaToH
if(WITH_PATOH)
  target_include_directories(xdecompose_mesh PRIVATE ${PATOH_INCLUDE_DIRS})
  target_link_libraries(xdecompose_mesh PRIVATE ${PATOH_LIBRARIES})
  target_compile_definitions(xdecompose_mesh PRIVATE USE_PATOH)
endif()

if(WITH_MPI)
  target_link_libraries(xdecompose_mesh PRIVATE MPI::MPI_Fortran)
endif()

if(WITH_HDF5)
  target_link_libraries(xdecompose_mesh PRIVATE ${HDF5_Fortran_LIBRARIES})
  target_include_directories(xdecompose_mesh PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
endif()

# Install
install(TARGETS xdecompose_mesh DESTINATION bin)

###############################################################################
# xdecompose_mesh_mpi (parallel version - only if WITH_MPI)
###############################################################################

if(WITH_MPI)
  add_executable(xdecompose_mesh_mpi
    $<TARGET_OBJECTS:decompose_mesh_par_obj>
    ${DECOMPOSE_MPI_MODULE_SOURCES}
    ${DECOMPOSE_MPI_SOURCES}
  )

  target_link_libraries(xdecompose_mesh_mpi PRIVATE
    specfem3d_shared
    MPI::MPI_Fortran
  )

  # Scotch
  if(WITH_SCOTCH)
    target_include_directories(xdecompose_mesh_mpi PRIVATE ${SCOTCH_INCLUDE_DIRS})
    target_link_libraries(xdecompose_mesh_mpi PRIVATE ${SCOTCH_LIBRARIES})
    target_compile_definitions(xdecompose_mesh_mpi PRIVATE USE_SCOTCH)
  endif()

  if(WITH_HDF5)
    target_link_libraries(xdecompose_mesh_mpi PRIVATE ${HDF5_Fortran_LIBRARIES})
    target_include_directories(xdecompose_mesh_mpi PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
  endif()

  # Install
  install(TARGETS xdecompose_mesh_mpi DESTINATION bin)
endif()
