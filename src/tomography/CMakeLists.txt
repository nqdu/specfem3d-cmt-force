# CMakeLists.txt for tomography tools
#
# Builds from src/tomography/:
#   - xadd_model_iso
#   - xmodel_update
#   - xsum_kernels
#   - xsum_preconditioned_kernels
#
# Builds from src/tomography/postprocess_sensitivity_kernels/:
#   - xclip_sem
#   - xcombine_sem
#   - xsmooth_sem
#   - xsmooth_sem_pde

###############################################################################
# Tomography module (shared by tomography executables)
###############################################################################

add_library(tomography_par_module OBJECT
  tomography_par.f90
)

target_link_libraries(tomography_par_module PUBLIC specfem3d_shared)

###############################################################################
# xadd_model_iso
###############################################################################

add_executable(xadd_model_iso
  add_model_iso.f90
  compute_kernel_integral.f90
  get_cg_direction.f90
  get_sd_direction.f90
  read_kernels.f90
  read_kernels_cg.f90
  read_model.f90
  read_parameters_tomo.f90
  write_gradient.f90
  write_new_model.f90
  write_new_model_perturbations.f90
)

target_link_libraries(xadd_model_iso PRIVATE
  tomography_par_module
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(WITH_MPI)
  target_link_libraries(xadd_model_iso PRIVATE MPI::MPI_Fortran)
endif()

if(HAS_GPU)
  target_link_libraries(xadd_model_iso PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_VTK)
  target_link_libraries(xadd_model_iso PRIVATE ${VTK_LIBRARIES})
endif()

###############################################################################
# xmodel_update
###############################################################################

add_executable(xmodel_update
  get_sd_direction.f90
  model_update.f90
  read_kernels.f90
  read_parameters_tomo.f90
  save_external_bin_m_up.f90
  write_gradient.f90
  write_new_model.f90
  write_new_model_perturbations.f90
)

target_link_libraries(xmodel_update PRIVATE
  tomography_par_module
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(WITH_MPI)
  target_link_libraries(xmodel_update PRIVATE MPI::MPI_Fortran)
endif()

if(HAS_GPU)
  target_link_libraries(xmodel_update PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_VTK)
  target_link_libraries(xmodel_update PRIVATE ${VTK_LIBRARIES})
endif()

###############################################################################
# xsum_kernels
###############################################################################

add_executable(xsum_kernels
  sum_kernels.f90
)

target_link_libraries(xsum_kernels PRIVATE
  tomography_par_module
  specfem3d_shared
)

if(WITH_MPI)
  target_link_libraries(xsum_kernels PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xsum_preconditioned_kernels
###############################################################################

add_executable(xsum_preconditioned_kernels
  sum_preconditioned_kernels.f90
)

target_link_libraries(xsum_preconditioned_kernels PRIVATE
  tomography_par_module
  specfem3d_shared
)

if(WITH_MPI)
  target_link_libraries(xsum_preconditioned_kernels PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# Postprocess module (shared by postprocess executables)
###############################################################################

add_library(postprocess_par_module OBJECT
  postprocess_sensitivity_kernels/postprocess_par.f90
)

target_link_libraries(postprocess_par_module PUBLIC specfem3d_shared)

###############################################################################
# xclip_sem
###############################################################################

add_executable(xclip_sem
  postprocess_sensitivity_kernels/parse_kernel_names.f90
  postprocess_sensitivity_kernels/clip_sem.f90
)

target_link_libraries(xclip_sem PRIVATE
  postprocess_par_module
  specfem3d_shared
)

if(WITH_MPI)
  target_link_libraries(xclip_sem PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xcombine_sem
###############################################################################

add_executable(xcombine_sem
  postprocess_sensitivity_kernels/parse_kernel_names.f90
  postprocess_sensitivity_kernels/combine_sem.f90
)

target_link_libraries(xcombine_sem PRIVATE
  postprocess_par_module
  specfem3d_shared
)

if(WITH_MPI)
  target_link_libraries(xcombine_sem PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xsmooth_sem
###############################################################################

add_executable(xsmooth_sem
  postprocess_sensitivity_kernels/parse_kernel_names.f90
  postprocess_sensitivity_kernels/smooth_sem.F90
)

target_link_libraries(xsmooth_sem PRIVATE
  postprocess_par_module
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(HAS_GPU)
  target_link_libraries(xsmooth_sem PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_MPI)
  target_link_libraries(xsmooth_sem PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xsmooth_sem_pde
###############################################################################

add_executable(xsmooth_sem_pde
  postprocess_sensitivity_kernels/parse_kernel_names.f90
  postprocess_sensitivity_kernels/smooth_sem_pde.F90
)

target_link_libraries(xsmooth_sem_pde PRIVATE
  postprocess_par_module
  specfem3d_solver_objects
  specfem3d_shared
  specfem3d_gpu
)

if(HAS_GPU)
  target_link_libraries(xsmooth_sem_pde PRIVATE ${GPU_LINK_LIBRARIES})
endif()

if(WITH_MPI)
  target_link_libraries(xsmooth_sem_pde PRIVATE MPI::MPI_Fortran)
endif()

if(WITH_ASDF)
  target_link_libraries(xsmooth_sem_pde PRIVATE ${ASDF_LIBRARIES})
endif()

###############################################################################
# Install
###############################################################################

install(TARGETS
  xadd_model_iso
  xmodel_update
  xsum_kernels
  xsum_preconditioned_kernels
  xclip_sem
  xcombine_sem
  xsmooth_sem
  xsmooth_sem_pde
  DESTINATION bin
)
