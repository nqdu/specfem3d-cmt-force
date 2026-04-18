###############################################################################
# generate_databases - builds xgenerate_databases executable
###############################################################################

# Module file (must be compiled first for module dependencies)
set(GENERATE_DATABASES_MODULE_SOURCES
  generate_databases_par.F90
)

# Main source files
set(GENERATE_DATABASES_SOURCES
  calc_jacobian.f90
  create_mass_matrices.f90
  create_regions_mesh.f90
  fault_generate_databases.f90
  finalize_databases.f90
  generate_databases.f90
  get_absorbing_boundary.F90
  get_coupling_surfaces.f90
  get_model.F90
  get_MPI.f90
  get_perm_color.f90
  lts_generate_databases.F90
  memory_eval.f90
  model_1d_cascadia.f90
  model_1d_prem.f90
  model_1d_socal.f90
  model_aniso.f90
  model_coupled.f90
  model_default.f90
  model_external_values.F90
  model_gll.f90
  model_ipati.f90
  model_salton_trough.f90
  model_scattering.F90
  model_tomography.f90
  model_vs30.f90
  pml_set_local_dampingcoeff.f90
  read_parameters.f90
  read_partition_files.f90
  save_arrays_solver.F90
  setup_color_perm.f90
  setup_mesh_adjacency.f90
  setup_mesh.f90
  wavefield_discontinuity_generate_databases.f90
)

# C source
set(GENERATE_DATABASES_C_SOURCES
  parse_sep.c
)

# HDF5 sources (always compiled, uses preprocessor for conditional HDF5 support)
set(GENERATE_DATABASES_HDF5_SOURCES
  read_partition_files_hdf5.F90
  save_arrays_solver_hdf5.F90
)

# ADIOS sources (conditional)
set(GENERATE_DATABASES_ADIOS_SOURCES
  read_partition_files_adios.F90
  save_arrays_solver_adios.F90
  save_moho_adios.F90
  model_gll_adios.F90
  model_ipati_adios.F90
)

# MPI-conditional sources
if(WITH_MPI)
  set(GENERATE_DATABASES_MPI_SOURCES model_sep.f90)
else()
  set(GENERATE_DATABASES_MPI_SOURCES model_sep_nompi.f90)
endif()

# Combine all sources
set(GENERATE_DATABASES_ALL_SOURCES
  ${GENERATE_DATABASES_MODULE_SOURCES}
  ${GENERATE_DATABASES_SOURCES}
  ${GENERATE_DATABASES_C_SOURCES}
  ${GENERATE_DATABASES_HDF5_SOURCES}
  ${GENERATE_DATABASES_MPI_SOURCES}
)

if(WITH_ADIOS OR WITH_ADIOS2)
  list(APPEND GENERATE_DATABASES_ALL_SOURCES ${GENERATE_DATABASES_ADIOS_SOURCES})
endif()

add_executable(xgenerate_databases ${GENERATE_DATABASES_ALL_SOURCES})

target_link_libraries(xgenerate_databases PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xgenerate_databases PRIVATE MPI::MPI_Fortran)
endif()

if(WITH_HDF5)
  target_link_libraries(xgenerate_databases PRIVATE ${HDF5_Fortran_LIBRARIES})
  target_include_directories(xgenerate_databases PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
endif()

if(WITH_ADIOS)
  target_link_libraries(xgenerate_databases PRIVATE ${ADIOS_LIBRARIES})
  target_include_directories(xgenerate_databases PRIVATE ${ADIOS_INCLUDE_DIRS})
elseif(WITH_ADIOS2)
  target_link_libraries(xgenerate_databases PRIVATE ${ADIOS2_LIBRARIES})
  target_include_directories(xgenerate_databases PRIVATE ${ADIOS2_INCLUDE_DIRS})
endif()

# Install
install(TARGETS xgenerate_databases DESTINATION bin)
