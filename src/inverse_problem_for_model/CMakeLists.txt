# CMakeLists.txt for inverse_problem_for_model
#
# Builds xinverse_problem_for_model executable (MPI only).
# Reuses specfem3D solver objects via specfem3d_solver_objects.

###############################################################################
# Source files from subdirectories
###############################################################################

set(INV_ADJOINT_SOURCE_SOURCES
  adjoint_source/signal_processing_mod.f90
  adjoint_source/rotations_mod.f90
  adjoint_source/adjoint_source_mod.f90
)

set(INV_INPUT_OUTPUT_SOURCES
  input_output/interpolation_mod.f90
  input_output/mesh_tools_mod.f90
  input_output/IO_model_mod.f90
  input_output/passive_imaging_format_mod.f90
  input_output/Teleseismic_IO.f90
  input_output/input_output_mod.f90
)

set(INV_INVERSION_SCHEME_SOURCES
  inversion_scheme/PrecondFWI_mod.f90
  inversion_scheme/iso_parameters.f90
  inversion_scheme/vti_parameters.f90
  inversion_scheme/family_parameter_mod.f90
  inversion_scheme/fwi_iteration_mod.f90
  inversion_scheme/inversion_scheme_mod.f90
)

set(INV_PARAMETERISATION_SOURCES
  parametrisation/elastic_isotropic_mod.f90
  parametrisation/parameters_subfamily_mod.f90
)

set(INV_PROJECTION_SOURCES
  projection_on_FD_grid/projection_on_FD_grid_mod.f90
)

set(INV_REGULARIZATION_SOURCES
  regularization/regularization_SEM_mod.f90
  regularization/regularization_FD_mod.f90
  regularization/regularization_interface.f90
)

set(INV_SPECFEM_INTERFACE_SOURCES
  specfem_interface/specfem_interface_mod.F90
)

set(INV_TOP_SOURCES
  inverse_problem_par.f90
  anisotropic_parametrisation_mod.f90
  elastic_tensor_tools_mod.f90
  inverse_problem_main.f90
  parallel_for_inverse_problem.f90
  program_inverse_problem.f90
)

###############################################################################
# Executable: xinverse_problem_for_model
###############################################################################

add_executable(xinverse_problem_for_model
  ${INV_TOP_SOURCES}
  ${INV_ADJOINT_SOURCE_SOURCES}
  ${INV_INPUT_OUTPUT_SOURCES}
  ${INV_INVERSION_SCHEME_SOURCES}
  ${INV_PARAMETERISATION_SOURCES}
  ${INV_PROJECTION_SOURCES}
  ${INV_REGULARIZATION_SOURCES}
  ${INV_SPECFEM_INTERFACE_SOURCES}
)

target_link_libraries(xinverse_problem_for_model PRIVATE
  specfem3d_solver_objects
  specfem3d_shared
)

if(WITH_MPI)
  target_link_libraries(xinverse_problem_for_model PRIVATE MPI::MPI_Fortran)
endif()

if(HAS_GPU)
  target_link_libraries(xinverse_problem_for_model PRIVATE specfem3d_gpu ${GPU_LINK_LIBRARIES})
endif()

if(WITH_VTK)
  target_link_libraries(xinverse_problem_for_model PRIVATE ${VTK_LIBRARIES})
endif()

if(WITH_ASDF)
  target_link_libraries(xinverse_problem_for_model PRIVATE ${ASDF_LIBRARIES})
endif()

###############################################################################
# Install
###############################################################################

install(TARGETS xinverse_problem_for_model DESTINATION bin)
