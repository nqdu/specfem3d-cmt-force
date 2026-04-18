# CMakeLists.txt for check_mesh_quality utilities
#
# Builds:
#   - xcheck_mesh_quality
#   - xconvert_skewness_to_angle

###############################################################################
# xcheck_mesh_quality
###############################################################################

add_executable(xcheck_mesh_quality
  check_mesh_quality.f90
)

target_link_libraries(xcheck_mesh_quality PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcheck_mesh_quality PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xconvert_skewness_to_angle
###############################################################################

add_executable(xconvert_skewness_to_angle
  convert_skewness_to_angle.f90
)

target_link_libraries(xconvert_skewness_to_angle PRIVATE specfem3d_shared)

###############################################################################
# Install
###############################################################################

install(TARGETS
  xcheck_mesh_quality
  xconvert_skewness_to_angle
  DESTINATION bin
)
