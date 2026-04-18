# CMakeLists.txt for auxiliary executables
#
# Builds:
#   - xcombine_surf_data
#   - xcombine_vol_data
#   - xcombine_vol_data_vtk
#   - xcombine_vol_data_vtu
#   - xcombine_vol_data_vtk_bin
#   - xconvolve_source_timefunction
#   - xdetect_duplicates_stations_file
#   - xcreate_movie_shakemap_AVS_DX_GMT
#   - xproject_and_combine_vol_data_on_regular_grid (MPI only)
#   - xcombine_vol_data_adios (ADIOS only)
#   - xcombine_vol_data_vtk_adios (ADIOS only)
#   - xcombine_vol_data_vtu_adios (ADIOS only)

###############################################################################
# xcombine_surf_data
###############################################################################

add_executable(xcombine_surf_data
  combine_surf_data.f90
)

target_link_libraries(xcombine_surf_data PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcombine_surf_data PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# combine_vol_data_impl module (shared by vol_data executables)
###############################################################################

add_library(combine_vol_data_impl_obj OBJECT
  combine_vol_data_impl.f90
)
target_link_libraries(combine_vol_data_impl_obj PUBLIC specfem3d_shared)

###############################################################################
# xcombine_vol_data
###############################################################################

add_executable(xcombine_vol_data
  combine_vol_data.F90
)

target_link_libraries(xcombine_vol_data PRIVATE combine_vol_data_impl_obj specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcombine_vol_data PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xcombine_vol_data_vtk
###############################################################################

add_executable(xcombine_vol_data_vtk
  combine_vol_data.F90
)

target_compile_definitions(xcombine_vol_data_vtk PRIVATE USE_VTK_INSTEAD_OF_MESH)
target_link_libraries(xcombine_vol_data_vtk PRIVATE combine_vol_data_impl_obj specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcombine_vol_data_vtk PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xcombine_vol_data_vtu
###############################################################################

add_executable(xcombine_vol_data_vtu
  combine_vol_data.F90
)

target_compile_definitions(xcombine_vol_data_vtu PRIVATE USE_VTU_INSTEAD_OF_MESH)
target_link_libraries(xcombine_vol_data_vtu PRIVATE combine_vol_data_impl_obj specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcombine_vol_data_vtu PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xcombine_vol_data_vtk_bin
###############################################################################

add_executable(xcombine_vol_data_vtk_bin
  combine_vol_data_vtk_binary.F90
  vtk_writer.c
)

target_link_libraries(xcombine_vol_data_vtk_bin PRIVATE combine_vol_data_impl_obj specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcombine_vol_data_vtk_bin PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xconvolve_source_timefunction
###############################################################################

add_executable(xconvolve_source_timefunction
  convolve_source_timefunction.f90
)

target_link_libraries(xconvolve_source_timefunction PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xconvolve_source_timefunction PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xdetect_duplicates_stations_file
###############################################################################

add_executable(xdetect_duplicates_stations_file
  detect_duplicates_stations_file.f90
)

target_link_libraries(xdetect_duplicates_stations_file PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xdetect_duplicates_stations_file PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################
# xcreate_movie_shakemap_AVS_DX_GMT
###############################################################################

add_executable(xcreate_movie_shakemap_AVS_DX_GMT
  create_movie_shakemap_AVS_DX_GMT.f90
)

target_link_libraries(xcreate_movie_shakemap_AVS_DX_GMT PRIVATE specfem3d_shared)

if(WITH_MPI)
  target_link_libraries(xcreate_movie_shakemap_AVS_DX_GMT PRIVATE MPI::MPI_Fortran)
endif()

###############################################################################


###############################################################################
# ADIOS variants (conditional)
###############################################################################

if(WITH_ADIOS OR WITH_ADIOS2)
  # xcombine_vol_data_adios
  add_executable(xcombine_vol_data_adios
    combine_vol_data.F90
    combine_vol_data_impl.f90
    combine_vol_data_adios_impl.f90
  )

  target_compile_definitions(xcombine_vol_data_adios PRIVATE USE_ADIOS_INSTEAD_OF_MESH)
  target_link_libraries(xcombine_vol_data_adios PRIVATE specfem3d_shared)

  if(WITH_MPI)
    target_link_libraries(xcombine_vol_data_adios PRIVATE MPI::MPI_Fortran)
  endif()
  if(WITH_ADIOS)
    target_link_libraries(xcombine_vol_data_adios PRIVATE ${ADIOS_LIBRARIES})
    target_include_directories(xcombine_vol_data_adios PRIVATE ${ADIOS_INCLUDE_DIRS})
  elseif(WITH_ADIOS2)
    target_link_libraries(xcombine_vol_data_adios PRIVATE ${ADIOS2_LIBRARIES})
    target_include_directories(xcombine_vol_data_adios PRIVATE ${ADIOS2_INCLUDE_DIRS})
  endif()

  # xcombine_vol_data_vtk_adios
  add_executable(xcombine_vol_data_vtk_adios
    combine_vol_data.F90
    combine_vol_data_impl.f90
    combine_vol_data_adios_impl.f90
  )

  target_compile_definitions(xcombine_vol_data_vtk_adios PRIVATE
    USE_ADIOS_INSTEAD_OF_MESH
    USE_VTK_INSTEAD_OF_MESH
  )
  target_link_libraries(xcombine_vol_data_vtk_adios PRIVATE specfem3d_shared)

  if(WITH_MPI)
    target_link_libraries(xcombine_vol_data_vtk_adios PRIVATE MPI::MPI_Fortran)
  endif()
  if(WITH_ADIOS)
    target_link_libraries(xcombine_vol_data_vtk_adios PRIVATE ${ADIOS_LIBRARIES})
    target_include_directories(xcombine_vol_data_vtk_adios PRIVATE ${ADIOS_INCLUDE_DIRS})
  elseif(WITH_ADIOS2)
    target_link_libraries(xcombine_vol_data_vtk_adios PRIVATE ${ADIOS2_LIBRARIES})
    target_include_directories(xcombine_vol_data_vtk_adios PRIVATE ${ADIOS2_INCLUDE_DIRS})
  endif()

  # xcombine_vol_data_vtu_adios
  add_executable(xcombine_vol_data_vtu_adios
    combine_vol_data.F90
    combine_vol_data_impl.f90
    combine_vol_data_adios_impl.f90
  )

  target_compile_definitions(xcombine_vol_data_vtu_adios PRIVATE
    USE_ADIOS_INSTEAD_OF_MESH
    USE_VTU_INSTEAD_OF_MESH
  )
  target_link_libraries(xcombine_vol_data_vtu_adios PRIVATE specfem3d_shared)

  if(WITH_MPI)
    target_link_libraries(xcombine_vol_data_vtu_adios PRIVATE MPI::MPI_Fortran)
  endif()
  if(WITH_ADIOS)
    target_link_libraries(xcombine_vol_data_vtu_adios PRIVATE ${ADIOS_LIBRARIES})
    target_include_directories(xcombine_vol_data_vtu_adios PRIVATE ${ADIOS_INCLUDE_DIRS})
  elseif(WITH_ADIOS2)
    target_link_libraries(xcombine_vol_data_vtu_adios PRIVATE ${ADIOS2_LIBRARIES})
    target_include_directories(xcombine_vol_data_vtu_adios PRIVATE ${ADIOS2_INCLUDE_DIRS})
  endif()

  install(TARGETS
    xcombine_vol_data_adios
    xcombine_vol_data_vtk_adios
    xcombine_vol_data_vtu_adios
    DESTINATION bin
  )
endif()

###############################################################################
# Install
###############################################################################

install(TARGETS
  xcombine_surf_data
  xcombine_vol_data
  xcombine_vol_data_vtk
  xcombine_vol_data_vtu
  xcombine_vol_data_vtk_bin
  xconvolve_source_timefunction
  xdetect_duplicates_stations_file
  xcreate_movie_shakemap_AVS_DX_GMT
  DESTINATION bin
)
