# CMakeLists.txt for specfem3D solver
#
# Builds:
#   1. specfem3d_solver_objects – OBJECT library with all solver objects
#      (reused by xstatic3D and xinverse_problem_for_model)
#   2. xspecfem3D – the main solver executable

###############################################################################
# Module sources (compiled first – other sources depend on these)
###############################################################################

set(SOLVER_MODULE_SOURCES
  specfem3D_par.F90
  asdf_data.f90
  pml_par.f90
)

###############################################################################
# Regular Fortran sources
###############################################################################

set(SOLVER_FORTRAN_SOURCES
  wavefield_discontinuity_solver.f90
  assemble_MPI_vector.f90
  check_stability.f90
  comp_source_time_function.f90
  compute_add_sources_acoustic.f90
  compute_add_sources_viscoelastic.F90
  compute_add_sources_poroelastic.f90
  compute_adj_source_frechet.f90
  compute_arrays_source.f90
  compute_boundary_kernel.f90
  compute_coupling_acoustic_el.f90
  compute_coupling_acoustic_po.f90
  compute_coupling_viscoelastic_ac.f90
  compute_coupling_viscoelastic_po.f90
  compute_coupling_poroelastic_ac.f90
  compute_coupling_poroelastic_el.f90
  compute_forces_acoustic_calling_routine.f90
  compute_forces_acoustic.F90
  compute_forces_viscoelastic_calling_routine.F90
  compute_forces_viscoelastic.F90
  compute_element_att_memory.f90
  compute_element_strain.F90
  compute_forces_poro_fluid_part.f90
  compute_forces_poroelastic_calling_routine.f90
  compute_forces_poro_solid_part.f90
  compute_gradient_in_acoustic.f90
  compute_interpolated_dva.f90
  compute_kernels.f90
  compute_seismograms.f90
  compute_stacey_acoustic.f90
  compute_stacey_viscoelastic.F90
  compute_stacey_poroelastic.f90
  compute_energy.f90
  convert_time.f90
  couple_with_injection.f90
  calendar.f90
  create_color_image.f90
  detect_mesh_surfaces.f90
  fault_solver_common.f90
  fault_solver_dynamic.f90
  fault_solver_kinematic.f90
  finalize_simulation.f90
  get_cmt.f90
  get_elevation.f90
  get_force.f90
  gravity_perturbation.f90
  initialize_simulation.F90
  iterate_time.F90
  iterate_time_undoatt.F90
  locate_MPI_slice.f90
  locate_point.f90
  locate_receivers.f90
  locate_source.F90
  lts_assemble_MPI_vector.F90
  lts_global_step.F90
  lts_iterate_time.F90
  lts_newmark_update.F90
  lts_setup.F90
  make_gravity.f90
  noise_tomography.f90
  pml_allocate_arrays.f90
  pml_output_VTKs.f90
  pml_compute_accel_contribution.f90
  pml_compute_memory_variables.f90
  prepare_attenuation.f90
  prepare_gpu.f90
  prepare_gravity.f90
  prepare_noise.f90
  prepare_optimized_arrays.F90
  prepare_timerun.F90
  prepare_wavefields.F90
  print_stf_file.f90
  read_external_stf.f90
  read_forward_arrays.f90
  read_mesh_databases.F90
  read_stations.f90
  save_adjoint_kernels.f90
  save_forward_arrays.f90
  setup_GLL_points.f90
  setup_movie_meshes.f90
  setup_sources_receivers.f90
  station_filter.f90
  surface_or_volume_integral.f90
  update_displacement_scheme.f90
  update_displacement_LDDRK.f90
  write_movie_output.F90
  write_output_ASCII_or_binary.f90
  write_output_SU.f90
  write_seismograms.f90
)

###############################################################################
# HDF5 sources (always compiled; use preprocessor guards internally)
###############################################################################

set(SOLVER_HDF5_SOURCES
  hdf5_io_server.F90
  read_mesh_databases_hdf5.F90
  write_movie_output_HDF5.F90
  write_output_HDF5.F90
)

###############################################################################
# Conditional: ADIOS / ADIOS2
###############################################################################

if(WITH_ADIOS OR WITH_ADIOS2)
  list(APPEND SOLVER_FORTRAN_SOURCES
    read_mesh_databases_adios.F90
    save_forward_arrays_adios.F90
    read_forward_arrays_adios.F90
    save_kernels_adios.F90
  )
endif()

###############################################################################
# Conditional: ASDF
###############################################################################

if(WITH_ASDF)
  list(APPEND SOLVER_FORTRAN_SOURCES
    write_output_ASDF.f90
    read_adjoint_sources_ASDF.f90
  )
endif()

###############################################################################
# Conditional: VTK
###############################################################################

if(WITH_VTK)
  set(SOLVER_VTK_SOURCES
    vtk_window.F90
    vtk_helper.cpp
  )
else()
  set(SOLVER_VTK_SOURCES
    vtk_window_stubs.c
  )
endif()

###############################################################################
# Object library: specfem3d_solver_objects
###############################################################################

add_library(specfem3d_solver_objects OBJECT
  ${SOLVER_MODULE_SOURCES}
  ${SOLVER_FORTRAN_SOURCES}
  ${SOLVER_HDF5_SOURCES}
  ${SOLVER_VTK_SOURCES}
)

target_link_libraries(specfem3d_solver_objects PUBLIC
  specfem3d_shared
  specfem3d_gpu
)

if(WITH_VTK)
  target_compile_definitions(specfem3d_solver_objects PUBLIC VTK_VIS)
  target_link_libraries(specfem3d_solver_objects PUBLIC ${VTK_LIBRARIES})
endif()

if(WITH_ADIOS)
  target_link_libraries(specfem3d_solver_objects PUBLIC ${ADIOS_LIBRARIES})
  target_include_directories(specfem3d_solver_objects PUBLIC ${ADIOS_INCLUDE_DIRS})
elseif(WITH_ADIOS2)
  target_link_libraries(specfem3d_solver_objects PUBLIC ${ADIOS2_LIBRARIES})
  target_include_directories(specfem3d_solver_objects PUBLIC ${ADIOS2_INCLUDE_DIRS})
endif()

if(WITH_ASDF)
  target_link_libraries(specfem3d_solver_objects PUBLIC ${ASDF_LIBRARIES})
endif()

if(WITH_HDF5)
  target_link_libraries(specfem3d_solver_objects PUBLIC ${HDF5_Fortran_LIBRARIES})
  target_include_directories(specfem3d_solver_objects PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
endif()

###############################################################################
# Executable: xspecfem3D
###############################################################################

add_executable(xspecfem3D specfem3D.F90)

target_link_libraries(xspecfem3D PRIVATE
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(WITH_MPI)
  target_link_libraries(xspecfem3D PRIVATE MPI::MPI_Fortran)
endif()

if(HAS_GPU)
  target_link_libraries(xspecfem3D PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_VTK)
  target_link_libraries(xspecfem3D PRIVATE ${VTK_LIBRARIES})
endif()

install(TARGETS xspecfem3D DESTINATION bin)
