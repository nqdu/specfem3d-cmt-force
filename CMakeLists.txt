#=====================================================================
#
#                          S p e c f e m 3 D
#                          -----------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                              CNRS, France
#                       and Princeton University, USA
#                 (there are currently many more authors!)
#                           (c) October 2017
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

cmake_minimum_required(VERSION 3.18)

project(specfem3d VERSION 4.1.1 LANGUAGES C CXX Fortran)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# CMake module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

###############################################################################
# Build Options
###############################################################################

option(WITH_MPI "Enable MPI support" ON)
option(DOUBLE_PRECISION "Use double precision" OFF)
option(ENABLE_DEBUG "Enable debug mode" OFF)
option(FORCE_VECTORIZATION "Force vectorization" OFF)
option(WITH_CUDA "Enable CUDA GPU support" OFF)
option(WITH_HIP "Enable HIP GPU support" OFF)
option(WITH_OPENMP "Enable OpenMP support" OFF)
option(WITH_VTK "Enable VTK visualization" OFF)
option(WITH_ADIOS "Enable ADIOS support" OFF)
option(WITH_ADIOS2 "Enable ADIOS2 support" OFF)
option(WITH_ASDF "Enable ASDF support" OFF)
option(WITH_HDF5 "Enable HDF5 support" OFF)
option(WITH_SCOTCH "Enable Scotch partitioning" ON)
option(USE_BUNDLED_SCOTCH "Use bundled Scotch library" ON)
option(WITH_METIS "Enable METIS partitioning" OFF)
option(WITH_PATOH "Enable PaToH partitioning" OFF)
option(WITH_PETSC "Enable PETSc linear solver backend" OFF)

set(CUDA_VERSION "" CACHE STRING "CUDA compute capability (e.g., cuda11)")
set(HIP_PLATFORM "" CACHE STRING "HIP platform (e.g., MI100, MI250)")
set(LOCAL_PATH_IS_ALSO_GLOBAL "false" CACHE STRING "Local path is also global")

###############################################################################
# Output Directories
###############################################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "Install prefix" FORCE)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

# Create output directories at configure time
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/obj)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/setup)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/OUTPUT_FILES)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/OUTPUT_FILES/DATABASES_MPI)

###############################################################################
# Compiler Detection & Flags
###############################################################################

# Include project compiler flags module (sets per-compiler Fortran/C flags)
include(cmake/CompilerFlags.cmake)

# Keep package-detection probes on a neutral Fortran toolchain. Project-specific
# flags are applied after MPI/OpenMP discovery so CMake's own test sources do not
# inherit strict diagnostics such as -fimplicit-none.
set(SPECFEM_Fortran_FLAGS_CHECK "${FLAGS_CHECK}")
set(SPECFEM_Fortran_FREE_FORM_FLAGS "")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(SPECFEM_Fortran_FREE_FORM_FLAGS "-ffree-form -ffree-line-length-none")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(SPECFEM_Fortran_FREE_FORM_FLAGS "-free")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(SPECFEM_Fortran_FREE_FORM_FLAGS "-Mfree")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(SPECFEM_Fortran_FREE_FORM_FLAGS "-f free")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  set(SPECFEM_Fortran_FREE_FORM_FLAGS "-ffree-form")
endif()

# Enable Fortran preprocessing globally
set(CMAKE_Fortran_PREPROCESS ON)

# Add include path for configured setup headers
include_directories(${CMAKE_BINARY_DIR}/setup)

###############################################################################
# Precision Settings
###############################################################################

if(DOUBLE_PRECISION)
  set(CUSTOM_REAL "SIZE_DOUBLE")
  set(CUSTOM_MPI_TYPE "MPI_DOUBLE_PRECISION")
else()
  set(CUSTOM_REAL "SIZE_REAL")
  set(CUSTOM_MPI_TYPE "MPI_REAL")
endif()

###############################################################################
# MPI
###############################################################################

if(WITH_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran C)
  message(STATUS "MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
  message(STATUS "MPI C compiler: ${MPI_C_COMPILER}")
  add_compile_definitions(WITH_MPI)
endif()

# Apply project-wide Fortran flags only after compiler-wrapper discovery.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SPECFEM_Fortran_FLAGS_CHECK} ${SPECFEM_Fortran_FREE_FORM_FLAGS}")

###############################################################################
# CUDA GPU Support
###############################################################################

set(HAS_GPU FALSE)
set(GPU_LINK_LIBRARIES "")

if(WITH_CUDA)
  enable_language(CUDA)

  # Map CUDA_VERSION string to SM architecture
  set(_cuda_arch "sm_60") # default fallback
  if(CUDA_VERSION STREQUAL "cuda5")
    set(_cuda_arch "sm_35")
  elseif(CUDA_VERSION STREQUAL "cuda6")
    set(_cuda_arch "sm_37")
  elseif(CUDA_VERSION STREQUAL "cuda7")
    set(_cuda_arch "sm_50")
  elseif(CUDA_VERSION STREQUAL "cuda8")
    set(_cuda_arch "sm_60")
  elseif(CUDA_VERSION STREQUAL "cuda9")
    set(_cuda_arch "sm_70")
  elseif(CUDA_VERSION STREQUAL "cuda10")
    set(_cuda_arch "sm_70")
  elseif(CUDA_VERSION STREQUAL "cuda11")
    set(_cuda_arch "sm_80")
  elseif(CUDA_VERSION STREQUAL "cuda12")
    set(_cuda_arch "sm_90")
  elseif(CUDA_VERSION STREQUAL "cuda13")
    set(_cuda_arch "sm_100")
  endif()

  # Extract numeric architecture from sm_XX
  string(REGEX REPLACE "sm_" "" _cuda_arch_num "${_cuda_arch}")
  set(CMAKE_CUDA_ARCHITECTURES ${_cuda_arch_num})

  message(STATUS "CUDA architecture: ${_cuda_arch} (CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES})")

  # Enable separable compilation for device linking
  set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

  set(HAS_GPU TRUE)
  add_compile_definitions(USE_CUDA)

  # Find CUDA toolkit for libraries
  find_package(CUDAToolkit REQUIRED)
  list(APPEND GPU_LINK_LIBRARIES CUDA::cudart)
endif()

###############################################################################
# HIP GPU Support
###############################################################################

if(WITH_HIP)
  find_package(hip REQUIRED)

  # Map HIP_PLATFORM string to GPU architecture
  set(_hip_arch "gfx906") # default fallback
  if(HIP_PLATFORM STREQUAL "MI8")
    set(_hip_arch "gfx803")
  elseif(HIP_PLATFORM STREQUAL "MI25")
    set(_hip_arch "gfx900")
  elseif(HIP_PLATFORM STREQUAL "MI50")
    set(_hip_arch "gfx906")
  elseif(HIP_PLATFORM STREQUAL "MI60")
    set(_hip_arch "gfx906")
  elseif(HIP_PLATFORM STREQUAL "MI100")
    set(_hip_arch "gfx908")
  elseif(HIP_PLATFORM STREQUAL "MI200")
    set(_hip_arch "gfx90a")
  elseif(HIP_PLATFORM STREQUAL "MI250")
    set(_hip_arch "gfx90a")
  elseif(HIP_PLATFORM STREQUAL "MI300")
    set(_hip_arch "gfx942")
  elseif(HIP_PLATFORM STREQUAL "MI350")
    set(_hip_arch "gfx950")
  endif()

  set(GPU_TARGETS "${_hip_arch}" CACHE STRING "HIP GPU architecture target")
  message(STATUS "HIP architecture: ${_hip_arch}")

  set(HAS_GPU TRUE)
  add_compile_definitions(USE_HIP)
  list(APPEND GPU_LINK_LIBRARIES hip::device)
endif()

###############################################################################
# OpenMP
###############################################################################

if(WITH_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran C)
  message(STATUS "OpenMP Fortran flags: ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

###############################################################################
# VTK
###############################################################################

if(WITH_VTK)
  find_package(VTK REQUIRED)
  message(STATUS "VTK version: ${VTK_VERSION}")
  set(HAVE_VTK TRUE)
else()
  set(HAVE_VTK FALSE)
endif()

###############################################################################
# ADIOS
###############################################################################

if(WITH_ADIOS)
  if(NOT WITH_MPI)
    message(FATAL_ERROR "Cannot compile with ADIOS without MPI. Enable WITH_MPI.")
  endif()
  find_package(ADIOS REQUIRED)
  message(STATUS "ADIOS found: ${ADIOS_LIBRARIES}")
endif()

###############################################################################
# ADIOS2
###############################################################################

if(WITH_ADIOS2)
  if(NOT WITH_MPI)
    message(FATAL_ERROR "Cannot compile with ADIOS2 without MPI. Enable WITH_MPI.")
  endif()
  find_package(ADIOS2 REQUIRED)
  message(STATUS "ADIOS2 found")
endif()

###############################################################################
# ASDF
###############################################################################

if(WITH_ASDF)
  find_package(ASDF REQUIRED)
  message(STATUS "ASDF found")
endif()

###############################################################################
# HDF5
###############################################################################

if(WITH_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS Fortran)
  message(STATUS "HDF5 Fortran libraries: ${HDF5_Fortran_LIBRARIES}")
endif()

###############################################################################
# Scotch Partitioning
###############################################################################

set(HAVE_SCOTCH FALSE)

if(WITH_SCOTCH)
  if(USE_BUNDLED_SCOTCH)
    # Build bundled Scotch from external_libs
    set(SCOTCH_SRC_DIR "${CMAKE_SOURCE_DIR}/external_libs/scotch_5.1.12b")
    set(SCOTCH_BUILD_DIR "${SCOTCH_SRC_DIR}/src")

    if(APPLE)
      set(SCOTCH_AR /usr/bin/ar)
      set(SCOTCH_RANLIB /usr/bin/ranlib)
    else()
      set(SCOTCH_AR ${CMAKE_AR})
      set(SCOTCH_RANLIB ${CMAKE_RANLIB})
    endif()

    message(STATUS "Using bundled Scotch from: ${SCOTCH_SRC_DIR}")

    if(APPLE)
      set(SCOTCH_MAKEFILE_INC "${SCOTCH_BUILD_DIR}/Make.inc/Makefile.inc.i686_mac_darwin8")
    else()
      set(SCOTCH_MAKEFILE_INC "${SCOTCH_BUILD_DIR}/Make.inc/Makefile.inc.x86-64_pc_linux2")
    endif()

    include(ExternalProject)
    ExternalProject_Add(scotch_bundled
      SOURCE_DIR "${SCOTCH_BUILD_DIR}"
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy
        "${SCOTCH_MAKEFILE_INC}"
        "${SCOTCH_BUILD_DIR}/Makefile.inc"
      BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -C "${SCOTCH_BUILD_DIR}"
        clean scotch
        "AR=${SCOTCH_AR}"
        "RANLIB=${SCOTCH_RANLIB}"
        "CCS=${CMAKE_C_COMPILER}"
        "CCD=${CMAKE_C_COMPILER}"
      BUILD_IN_SOURCE TRUE
      INSTALL_COMMAND ""
      BUILD_BYPRODUCTS
        "${SCOTCH_SRC_DIR}/lib/libscotch.a"
        "${SCOTCH_SRC_DIR}/lib/libscotcherr.a"
        "${SCOTCH_SRC_DIR}/lib/libscotcherrexit.a"
    )

    # Create imported targets for bundled scotch
    add_library(scotch_lib STATIC IMPORTED)
    set_target_properties(scotch_lib PROPERTIES
      IMPORTED_LOCATION "${SCOTCH_SRC_DIR}/lib/libscotch.a"
    )
    add_dependencies(scotch_lib scotch_bundled)

    add_library(scotcherr_lib STATIC IMPORTED)
    set_target_properties(scotcherr_lib PROPERTIES
      IMPORTED_LOCATION "${SCOTCH_SRC_DIR}/lib/libscotcherr.a"
    )
    add_dependencies(scotcherr_lib scotch_bundled)

    set(SCOTCH_LIBRARIES scotch_lib scotcherr_lib)
    set(SCOTCH_INCLUDE_DIRS "${SCOTCH_SRC_DIR}/include")
    set(HAVE_SCOTCH TRUE)

  else()
    # Find system-installed Scotch
    find_package(Scotch REQUIRED)
    message(STATUS "Scotch libraries: ${SCOTCH_LIBRARIES}")
    message(STATUS "Scotch include dirs: ${SCOTCH_INCLUDE_DIRS}")
    set(HAVE_SCOTCH TRUE)
  endif()
endif()

###############################################################################
# METIS Partitioning
###############################################################################

if(WITH_METIS)
  find_package(METIS REQUIRED)
  message(STATUS "METIS found")
endif()

###############################################################################
# PaToH Partitioning
###############################################################################

if(WITH_PATOH)
  find_package(PaToH REQUIRED)
  message(STATUS "PaToH found")
endif()

###############################################################################
# PETSc Linear Solver
###############################################################################

if(WITH_PETSC)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc petsc)
  message(STATUS "PETSc version: ${PETSC_VERSION}")
  message(STATUS "PETSc include dirs: ${PETSC_INCLUDE_DIRS}")
  # Create an alias so sub-directories can use PETSc::PETSc
  if(NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc ALIAS PkgConfig::PETSC)
  endif()
  add_compile_definitions(WITH_PETSC)
endif()

###############################################################################
# Force Vectorization
###############################################################################

if(FORCE_VECTORIZATION)
  add_compile_definitions(FORCE_VECTORIZATION)
endif()

###############################################################################
# Debug Mode
###############################################################################

if(ENABLE_DEBUG)
  add_compile_definitions(DEBUG)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -fcheck=all -fbacktrace -Wall -Wextra")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -check all -traceback -warn all")
  endif()
endif()

###############################################################################
# Header Checks for config.h
###############################################################################

include(CheckIncludeFile)

check_include_file("emmintrin.h" HAVE_EMMINTRIN)
check_include_file("xmmintrin.h" HAVE_XMMINTRIN)
check_include_file("err.h" HAVE_ERR)
check_include_file("inttypes.h" HAVE_INTTYPES_H)
check_include_file("stdint.h" HAVE_STDINT_H)
check_include_file("stdio.h" HAVE_STDIO_H)
check_include_file("stdlib.h" HAVE_STDLIB_H)
check_include_file("string.h" HAVE_STRING_H)
check_include_file("strings.h" HAVE_STRINGS_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("sys/types.h" HAVE_SYS_TYPES_H)
check_include_file("unistd.h" HAVE_UNISTD_H)

# Check for pthreads
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT)
  set(HAVE_PTHREAD TRUE)
endif()

# STDC_HEADERS: assume true on modern systems
set(STDC_HEADERS TRUE)

###############################################################################
# Git Version Information
###############################################################################

execute_process(
  COMMAND git describe --tags --dirty --always
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_PACKAGE_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
execute_process(
  COMMAND git log -1 --format=%H
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
execute_process(
  COMMAND git log -1 --format=%ci
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_DATE_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

if(NOT GIT_PACKAGE_VERSION)
  set(GIT_PACKAGE_VERSION "${PROJECT_VERSION}")
endif()
if(NOT GIT_COMMIT_VERSION)
  set(GIT_COMMIT_VERSION "unknown")
endif()
if(NOT GIT_DATE_VERSION)
  set(GIT_DATE_VERSION "unknown")
endif()

message(STATUS "Git version: ${GIT_PACKAGE_VERSION}")
message(STATUS "Git commit:  ${GIT_COMMIT_VERSION}")
message(STATUS "Git date:    ${GIT_DATE_VERSION}")

# Set variables for config.h.cmake.in
set(SPECFEM3D_GIT_REVISION "${GIT_PACKAGE_VERSION}")
set(SPECFEM3D_GIT_HASH "${GIT_COMMIT_VERSION}")
set(SPECFEM3D_GIT_DATE "${GIT_DATE_VERSION}")

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE SPECFEM3D_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT SPECFEM3D_GIT_BRANCH)
  set(SPECFEM3D_GIT_BRANCH "unknown")
endif()

###############################################################################
# Configure Headers
###############################################################################

# Detect Fortran name-mangling convention and define FC_FUNC / FC_FUNC_
# so that C sources (e.g. force_ftz.c) can call Fortran routines portably.
# FortranCInterface_HEADER with MACRO_NAMESPACE "FC_" produces:
#   FC_GLOBAL(name,NAME)  -- for identifiers WITHOUT underscores
#   FC_GLOBAL_(name,NAME) -- for identifiers WITH    underscores
# config.h.cmake.in then aliases these to the expected FC_FUNC / FC_FUNC_.
include(FortranCInterface)
FortranCInterface_HEADER(
  ${CMAKE_BINARY_DIR}/setup/FC_mangling.h
  MACRO_NAMESPACE "FC_"
)

# Build the configure flags string
set(CONFIGURE_FLAGS "FC=${CMAKE_Fortran_COMPILER} FCFLAGS=${CMAKE_Fortran_FLAGS}")
set(PACKAGE_STRING "SPECFEM3D Cartesian ${PROJECT_VERSION}")

# setup/constants.h.in -> ${CMAKE_BINARY_DIR}/setup/constants.h
# Uses @CUSTOM_REAL@ and @LOCAL_PATH_IS_ALSO_GLOBAL@
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/constants.h.in
  ${CMAKE_BINARY_DIR}/setup/constants.h
  @ONLY
)

# setup/precision.h.in -> ${CMAKE_BINARY_DIR}/setup/precision.h
# Uses @CUSTOM_MPI_TYPE@
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/precision.h.in
  ${CMAKE_BINARY_DIR}/setup/precision.h
  @ONLY
)

# setup/config.fh.in -> ${CMAKE_BINARY_DIR}/setup/config.fh
# Uses @CONFIGURE_FLAGS@ and @PACKAGE_STRING@
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/config.fh.in
  ${CMAKE_BINARY_DIR}/setup/config.fh
  @ONLY
)

# setup/constants_tomography.h.in -> ${CMAKE_BINARY_DIR}/setup/constants_tomography.h
# Uses @configure_input@ which expands but is inside a comment so harmless
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/constants_tomography.h.in
  ${CMAKE_BINARY_DIR}/setup/constants_tomography.h
  @ONLY
)

# setup/version.fh.in -> ${CMAKE_BINARY_DIR}/setup/version.fh
# Uses @GIT_PACKAGE_VERSION@, @GIT_COMMIT_VERSION@, @GIT_DATE_VERSION@
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/version.fh.in
  ${CMAKE_BINARY_DIR}/setup/version.fh
  @ONLY
)

# setup/config.h.cmake.in -> ${CMAKE_BINARY_DIR}/setup/config.h
# Uses #cmakedefine for feature detection macros
configure_file(
  ${CMAKE_SOURCE_DIR}/setup/config.h.cmake.in
  ${CMAKE_BINARY_DIR}/setup/config.h
)

###############################################################################
# Source Subdirectories
###############################################################################

add_subdirectory(src/shared)
add_subdirectory(src/gpu)
add_subdirectory(src/meshfem3D)
add_subdirectory(src/generate_databases)
add_subdirectory(src/specfem3D)
add_subdirectory(src/static)
add_subdirectory(src/decompose_mesh)
add_subdirectory(src/auxiliaries)
add_subdirectory(src/check_mesh_quality)
add_subdirectory(src/tomography)



###############################################################################
# Configuration Summary
###############################################################################

message(STATUS "")
message(STATUS "=============================================")
message(STATUS "  SPECFEM3D Cartesian ${PROJECT_VERSION} - Configuration Summary")
message(STATUS "=============================================")
message(STATUS "")
message(STATUS "  Build type:           ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Fortran compiler:     ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Fortran compiler ID:  ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  C compiler:           ${CMAKE_C_COMPILER}")
message(STATUS "  CXX compiler:         ${CMAKE_CXX_COMPILER}")
message(STATUS "")
message(STATUS "  Precision:")
if(DOUBLE_PRECISION)
  message(STATUS "    Double precision:   YES")
else()
  message(STATUS "    Double precision:   NO (single precision)")
endif()
message(STATUS "")
message(STATUS "  Parallelism:")
if(WITH_MPI)
  message(STATUS "    MPI:                YES")
  message(STATUS "      MPI Fortran:      ${MPI_Fortran_COMPILER}")
else()
  message(STATUS "    MPI:                NO")
endif()
if(WITH_OPENMP)
  message(STATUS "    OpenMP:             YES")
else()
  message(STATUS "    OpenMP:             NO")
endif()
message(STATUS "")
message(STATUS "  GPU:")
if(WITH_CUDA)
  message(STATUS "    CUDA:               YES (${CUDA_VERSION}, arch=${_cuda_arch})")
else()
  message(STATUS "    CUDA:               NO")
endif()
if(WITH_HIP)
  message(STATUS "    HIP:                YES (${HIP_PLATFORM}, arch=${_hip_arch})")
else()
  message(STATUS "    HIP:                NO")
endif()
message(STATUS "")
message(STATUS "  Partitioning:")
if(WITH_SCOTCH)
  if(USE_BUNDLED_SCOTCH)
    message(STATUS "    Scotch:             YES (bundled)")
  else()
    message(STATUS "    Scotch:             YES (system)")
  endif()
else()
  message(STATUS "    Scotch:             NO")
endif()
if(WITH_METIS)
  message(STATUS "    METIS:              YES")
else()
  message(STATUS "    METIS:              NO")
endif()
if(WITH_PATOH)
  message(STATUS "    PaToH:              YES")
else()
  message(STATUS "    PaToH:              NO")
endif()
message(STATUS "")
message(STATUS "  I/O Libraries:")
if(WITH_ADIOS)
  message(STATUS "    ADIOS:              YES")
else()
  message(STATUS "    ADIOS:              NO")
endif()
if(WITH_ADIOS2)
  message(STATUS "    ADIOS2:             YES")
else()
  message(STATUS "    ADIOS2:             NO")
endif()
if(WITH_ASDF)
  message(STATUS "    ASDF:               YES")
else()
  message(STATUS "    ASDF:               NO")
endif()
if(WITH_HDF5)
  message(STATUS "    HDF5:               YES")
else()
  message(STATUS "    HDF5:               NO")
endif()
message(STATUS "")
message(STATUS "  Other:")
if(WITH_VTK)
  message(STATUS "    VTK:                YES (${VTK_VERSION})")
else()
  message(STATUS "    VTK:                NO")
endif()
if(FORCE_VECTORIZATION)
  message(STATUS "    Vectorization:      FORCED")
else()
  message(STATUS "    Vectorization:      NO")
endif()
if(ENABLE_DEBUG)
  message(STATUS "    Debug mode:         YES")
else()
  message(STATUS "    Debug mode:         NO")
endif()
message(STATUS "    Local=Global path:  ${LOCAL_PATH_IS_ALSO_GLOBAL}")
message(STATUS "")
message(STATUS "  Git:")
message(STATUS "    Version:            ${GIT_PACKAGE_VERSION}")
message(STATUS "    Commit:             ${GIT_COMMIT_VERSION}")
message(STATUS "    Date:               ${GIT_DATE_VERSION}")
message(STATUS "    Branch:             ${SPECFEM3D_GIT_BRANCH}")
message(STATUS "")
message(STATUS "=============================================")
message(STATUS "")
