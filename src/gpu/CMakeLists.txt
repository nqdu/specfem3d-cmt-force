# GPU library - conditional on CUDA or HIP
#
# Builds OBJECT library specfem3d_gpu with GPU acceleration code.
# When neither CUDA nor HIP is enabled, compiles a stub library instead.

###############################################################################
# Source file lists
###############################################################################

# Main GPU source files (.cu in src/gpu/)
set(GPU_MAIN_SOURCES
  assemble_MPI_scalar_cuda.cu
  assemble_MPI_vector_cuda.cu
  check_fields_cuda.cu
  compute_add_sources_acoustic_cuda.cu
  compute_add_sources_viscoelastic_cuda.cu
  compute_coupling_cuda.cu
  compute_forces_acoustic_cuda.cu
  compute_forces_viscoelastic_cuda.cu
  compute_kernels_cuda.cu
  compute_stacey_acoustic_cuda.cu
  compute_stacey_viscoelastic_cuda.cu
  compute_strain_cuda.cu
  fault_solver_dynamics.cu
  helper_functions.cu
  initialize_gpu.cu
  lts_assembly_mpi_cuda.cu
  lts_compute_forces_viscoelastic_cuda.cu
  noise_tomography_cuda.cu
  pml_compute_accel_cuda.cu
  prepare_mesh_constants_cuda.cu
  smooth_cuda.cu
  smooth_pde_cuda.cu
  transfer_fields_cuda.cu
  update_displacement_cuda.cu
  wavefield_discontinuity_cuda.cu
  write_seismograms_cuda.cu
)

# GPU C sources in src/gpu/ (also built in Make via gpu_specfem3D_OBJECTS)
set(GPU_MAIN_C_SOURCES
  save_and_compare_cpu_vs_gpu.c
)

# CUDA kernel sources (.cu in src/gpu/kernels/)
set(KERNEL_CUDA_SOURCES
  kernels/Kernel_2_acoustic_impl.cu
  kernels/Kernel_2_viscoelastic_impl.cu
  kernels/UpdateDispVeloc_kernel.cu
  kernels/UpdatePotential_kernel.cu
  kernels/add_source_main_rec_noise_cuda_kernel.cu
  kernels/add_sources_ac_SIM_TYPE_2_OR_3_kernel.cu
  kernels/add_sources_el_SIM_TYPE_2_OR_3_kernel.cu
  kernels/assemble_boundary_accel_on_device.cu
  kernels/assemble_boundary_potential_on_device.cu
  kernels/compute_acoustic_seismogram_kernel.cu
  kernels/compute_add_sources_acoustic_kernel.cu
  kernels/compute_add_sources_kernel.cu
  kernels/compute_coupling_acoustic_el_kernel.cu
  kernels/compute_coupling_elastic_ac_kernel.cu
  kernels/compute_coupling_ocean_cuda_kernel.cu
  kernels/compute_dynamic_fault_cuda.cu
  kernels/compute_element_strain_cudakernel.cu
  kernels/compute_elastic_seismogram_kernel.cu
  kernels/compute_kernels_acoustic_kernel.cu
  kernels/compute_kernels_ani_cudakernel.cu
  kernels/compute_kernels_cudakernel.cu
  kernels/compute_kernels_hess_ac_cudakernel.cu
  kernels/compute_kernels_hess_el_cudakernel.cu
  kernels/compute_kernels_strength_noise_cuda_kernel.cu
  kernels/compute_stacey_acoustic_kernel.cu
  kernels/compute_stacey_elastic_kernel.cu
  kernels/enforce_free_surface_cuda_kernel.cu
  kernels/get_maximum_field_kernel.cu
  kernels/get_maximum_vector_kernel.cu
  kernels/kernel_3_accel_cuda_device.cu
  kernels/kernel_3_acoustic_cuda_device.cu
  kernels/kernel_3_cuda_device.cu
  kernels/kernel_3_veloc_cuda_device.cu
  kernels/lts_assembly_mpi_kernel.cu
  kernels/lts_compute_forces_viscoelastic_kernel.cu
  kernels/noise_read_add_surface_movie_cuda_kernel.cu
  kernels/pml_impose_boundary_condition_cuda_kernel.cu
  kernels/pml_kernel_2_viscoelastic_impl.cu
  kernels/prepare_boundary_accel_on_device.cu
  kernels/prepare_boundary_potential_on_device.cu
  kernels/process_smooth.cu
  kernels/smooth_pde_kernels.cu
  kernels/synchronize_boundary_accel_on_device.cu
  kernels/transfer_surface_to_host_kernel.cu
  kernels/wavefield_discontinuity_kernel.cu
)

# HIP kernel sources (.cpp in src/gpu/kernels/)
set(KERNEL_HIP_SOURCES
  kernels/Kernel_2_acoustic_impl.cpp
  kernels/Kernel_2_viscoelastic_impl.cpp
  kernels/UpdateDispVeloc_kernel.cpp
  kernels/UpdatePotential_kernel.cpp
  kernels/add_source_main_rec_noise_cuda_kernel.cpp
  kernels/add_sources_ac_SIM_TYPE_2_OR_3_kernel.cpp
  kernels/add_sources_el_SIM_TYPE_2_OR_3_kernel.cpp
  kernels/assemble_boundary_accel_on_device.cpp
  kernels/assemble_boundary_potential_on_device.cpp
  kernels/compute_acoustic_seismogram_kernel.cpp
  kernels/compute_add_sources_acoustic_kernel.cpp
  kernels/compute_add_sources_kernel.cpp
  kernels/compute_coupling_acoustic_el_kernel.cpp
  kernels/compute_coupling_elastic_ac_kernel.cpp
  kernels/compute_coupling_ocean_cuda_kernel.cpp
  kernels/compute_dynamic_fault_cuda.cpp
  kernels/compute_element_strain_cudakernel.cpp
  kernels/compute_elastic_seismogram_kernel.cpp
  kernels/compute_kernels_acoustic_kernel.cpp
  kernels/compute_kernels_ani_cudakernel.cpp
  kernels/compute_kernels_cudakernel.cpp
  kernels/compute_kernels_hess_ac_cudakernel.cpp
  kernels/compute_kernels_hess_el_cudakernel.cpp
  kernels/compute_kernels_strength_noise_cuda_kernel.cpp
  kernels/compute_stacey_acoustic_kernel.cpp
  kernels/compute_stacey_elastic_kernel.cpp
  kernels/enforce_free_surface_cuda_kernel.cpp
  kernels/get_maximum_field_kernel.cpp
  kernels/get_maximum_vector_kernel.cpp
  kernels/kernel_3_accel_cuda_device.cpp
  kernels/kernel_3_acoustic_cuda_device.cpp
  kernels/kernel_3_cuda_device.cpp
  kernels/kernel_3_veloc_cuda_device.cpp
  kernels/lts_assembly_mpi_kernel.cpp
  kernels/lts_compute_forces_viscoelastic_kernel.cpp
  kernels/noise_read_add_surface_movie_cuda_kernel.cpp
  kernels/pml_impose_boundary_condition_cuda_kernel.cpp
  kernels/pml_kernel_2_viscoelastic_impl.cpp
  kernels/prepare_boundary_accel_on_device.cpp
  kernels/prepare_boundary_potential_on_device.cpp
  kernels/process_smooth.cpp
  kernels/smooth_pde_kernels.cpp
  kernels/synchronize_boundary_accel_on_device.cpp
  kernels/transfer_surface_to_host_kernel.cpp
  kernels/wavefield_discontinuity_kernel.cpp
)

###############################################################################
# Build targets
###############################################################################

if(WITH_CUDA)
  # -------------------------------------------------------------------
  # CUDA build: .cu files compiled by nvcc via CMake's CUDA language support
  # -------------------------------------------------------------------
  add_library(specfem3d_gpu OBJECT
    ${GPU_MAIN_SOURCES}
    ${GPU_MAIN_C_SOURCES}
    ${KERNEL_CUDA_SOURCES}
  )

  target_include_directories(specfem3d_gpu PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_BINARY_DIR}/setup
  )

  target_compile_definitions(specfem3d_gpu PUBLIC USE_CUDA)

  if(WITH_MPI)
    target_compile_definitions(specfem3d_gpu PUBLIC WITH_MPI)
    target_include_directories(specfem3d_gpu PUBLIC ${MPI_C_INCLUDE_DIRS})
  endif()

  # Kernel files need mesh_constants_gpu.h force-included so device code
  # picks up the constant declarations.
  set_source_files_properties(
    ${KERNEL_CUDA_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-include;${CMAKE_CURRENT_SOURCE_DIR}/mesh_constants_gpu.h"
  )

elseif(WITH_HIP)
  # -------------------------------------------------------------------
  # HIP build: main .cu files compiled by hipcc (language set to HIP),
  # kernel .cpp wrappers compiled normally by hipcc.
  # -------------------------------------------------------------------
  add_library(specfem3d_gpu OBJECT
    ${GPU_MAIN_SOURCES}
    ${GPU_MAIN_C_SOURCES}
    ${KERNEL_HIP_SOURCES}
  )

  target_include_directories(specfem3d_gpu PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels
    ${CMAKE_BINARY_DIR}/setup
  )

  target_compile_definitions(specfem3d_gpu PUBLIC USE_HIP)
  target_link_libraries(specfem3d_gpu PUBLIC hip::device)

  if(WITH_MPI)
    target_compile_definitions(specfem3d_gpu PUBLIC WITH_MPI)
    target_include_directories(specfem3d_gpu PUBLIC ${MPI_C_INCLUDE_DIRS})
  endif()

  # Tell CMake to treat the .cu main sources as HIP language so hipcc
  # compiles them instead of trying to use a C++ compiler.
  set_source_files_properties(
    ${GPU_MAIN_SOURCES}
    PROPERTIES LANGUAGE HIP
  )

  # Kernel .cpp files also need mesh_constants_gpu.h force-included.
  set_source_files_properties(
    ${KERNEL_HIP_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-include;${CMAKE_CURRENT_SOURCE_DIR}/mesh_constants_gpu.h"
  )

else()
  # -------------------------------------------------------------------
  # No GPU: compile C stubs that provide empty implementations
  # -------------------------------------------------------------------
  add_library(specfem3d_gpu OBJECT
    specfem3D_gpu_cuda_method_stubs.c
  )

  target_include_directories(specfem3d_gpu PUBLIC
    ${CMAKE_BINARY_DIR}/setup
  )

  if(WITH_MPI)
    target_include_directories(specfem3d_gpu PUBLIC ${MPI_C_INCLUDE_DIRS})
  endif()
endif()
