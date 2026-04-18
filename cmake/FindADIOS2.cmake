# FindADIOS2.cmake
# -----------------
# Find the ADIOS2 I/O library (version 2.x).
#
# This module first attempts to find ADIOS2 via its CMake config package
# (the preferred method). If that fails, it falls back to pkg-config, and
# finally to manual header/library search.
#
# Hints:
#   ADIOS2_DIR       - Root directory of ADIOS2 installation (env or CMake var)
#   ADIOS2_ROOT      - Alternative root directory hint
#
# Result variables:
#   ADIOS2_FOUND         - True if ADIOS2 was found
#   ADIOS2_INCLUDE_DIRS  - Include directories for ADIOS2
#   ADIOS2_LIBRARIES     - Libraries to link against
#   ADIOS2_VERSION       - ADIOS2 version string (if available)

include(FindPackageHandleStandardArgs)

set(_adios2_FOUND FALSE)

# Collect search hints
set(_adios2_SEARCH_DIRS "")
if(ADIOS2_DIR)
    list(APPEND _adios2_SEARCH_DIRS "${ADIOS2_DIR}")
endif()
if(DEFINED ENV{ADIOS2_DIR})
    list(APPEND _adios2_SEARCH_DIRS "$ENV{ADIOS2_DIR}")
endif()
if(ADIOS2_ROOT)
    list(APPEND _adios2_SEARCH_DIRS "${ADIOS2_ROOT}")
endif()
if(DEFINED ENV{ADIOS2_ROOT})
    list(APPEND _adios2_SEARCH_DIRS "$ENV{ADIOS2_ROOT}")
endif()

# --- Strategy 1: CMake config mode ---
# ADIOS2 ships its own ADIOS2Config.cmake
set(_adios2_config_search_paths "")
foreach(_dir ${_adios2_SEARCH_DIRS})
    list(APPEND _adios2_config_search_paths
        "${_dir}"
        "${_dir}/lib/cmake/adios2"
        "${_dir}/lib64/cmake/adios2"
    )
endforeach()

find_package(ADIOS2 CONFIG QUIET
    HINTS ${_adios2_config_search_paths}
)

if(ADIOS2_FOUND)
    # CMake config mode found it. Extract information if targets exist.
    set(_adios2_FOUND TRUE)

    # ADIOS2 provides imported targets like adios2::adios2_f or adios2::cxx11_mpi
    # Collect what is available
    set(ADIOS2_LIBRARIES "")
    set(ADIOS2_INCLUDE_DIRS "")

    # Prefer the Fortran MPI target for SPECFEM3D
    foreach(_target
        adios2::adios2_f_mpi
        adios2::adios2_f
        adios2::cxx11_mpi
        adios2::cxx11
        adios2::adios2
    )
        if(TARGET ${_target})
            list(APPEND ADIOS2_LIBRARIES "${_target}")
            get_target_property(_inc ${_target} INTERFACE_INCLUDE_DIRECTORIES)
            if(_inc)
                list(APPEND ADIOS2_INCLUDE_DIRS ${_inc})
            endif()
            break()
        endif()
    endforeach()

    if(ADIOS2_INCLUDE_DIRS)
        list(REMOVE_DUPLICATES ADIOS2_INCLUDE_DIRS)
    endif()

    if(ADIOS2_VERSION)
        # Already set by the config package
    endif()
endif()

# --- Strategy 2: pkg-config ---
if(NOT _adios2_FOUND)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        # Set PKG_CONFIG_PATH if we have hints
        set(_saved_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
        foreach(_dir ${_adios2_SEARCH_DIRS})
            set(ENV{PKG_CONFIG_PATH} "${_dir}/lib/pkgconfig:${_dir}/lib64/pkgconfig:$ENV{PKG_CONFIG_PATH}")
        endforeach()

        pkg_check_modules(_adios2_pkg QUIET adios2)

        # Restore PKG_CONFIG_PATH
        set(ENV{PKG_CONFIG_PATH} "${_saved_pkg_config_path}")

        if(_adios2_pkg_FOUND)
            set(_adios2_FOUND TRUE)
            set(ADIOS2_INCLUDE_DIRS "${_adios2_pkg_INCLUDE_DIRS}")
            set(ADIOS2_LIBRARIES "${_adios2_pkg_LIBRARIES}")
            set(ADIOS2_VERSION "${_adios2_pkg_VERSION}")

            # Resolve library names to full paths if possible
            set(_adios2_resolved_libs "")
            foreach(_lib ${_adios2_pkg_LIBRARIES})
                find_library(_adios2_lib_${_lib}
                    NAMES ${_lib}
                    HINTS ${_adios2_pkg_LIBRARY_DIRS}
                )
                if(_adios2_lib_${_lib})
                    list(APPEND _adios2_resolved_libs "${_adios2_lib_${_lib}}")
                else()
                    list(APPEND _adios2_resolved_libs "${_lib}")
                endif()
            endforeach()
            set(ADIOS2_LIBRARIES "${_adios2_resolved_libs}")
        endif()
    endif()
endif()

# --- Strategy 3: Manual search (last resort) ---
if(NOT _adios2_FOUND)
    find_path(ADIOS2_INCLUDE_DIR
        NAMES adios2.h
        HINTS
            ${_adios2_SEARCH_DIRS}
        PATH_SUFFIXES
            include
        DOC "Directory containing adios2.h"
    )

    find_library(ADIOS2_CORE_LIBRARY
        NAMES adios2_core adios2
        HINTS
            ${_adios2_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS2 core library"
    )

    find_library(ADIOS2_CORE_MPI_LIBRARY
        NAMES adios2_core_mpi
        HINTS
            ${_adios2_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS2 core MPI library"
    )

    # Fortran bindings
    find_library(ADIOS2_FORTRAN_LIBRARY
        NAMES adios2_fortran adios2_f
        HINTS
            ${_adios2_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS2 Fortran library"
    )

    find_library(ADIOS2_FORTRAN_MPI_LIBRARY
        NAMES adios2_fortran_mpi adios2_f_mpi
        HINTS
            ${_adios2_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS2 Fortran MPI library"
    )

    set(ADIOS2_INCLUDE_DIRS "")
    if(ADIOS2_INCLUDE_DIR)
        set(ADIOS2_INCLUDE_DIRS "${ADIOS2_INCLUDE_DIR}")
    endif()

    # Collect libraries in link order
    set(ADIOS2_LIBRARIES "")
    if(ADIOS2_FORTRAN_MPI_LIBRARY)
        list(APPEND ADIOS2_LIBRARIES "${ADIOS2_FORTRAN_MPI_LIBRARY}")
    elseif(ADIOS2_FORTRAN_LIBRARY)
        list(APPEND ADIOS2_LIBRARIES "${ADIOS2_FORTRAN_LIBRARY}")
    endif()
    if(ADIOS2_CORE_MPI_LIBRARY)
        list(APPEND ADIOS2_LIBRARIES "${ADIOS2_CORE_MPI_LIBRARY}")
    endif()
    if(ADIOS2_CORE_LIBRARY)
        list(APPEND ADIOS2_LIBRARIES "${ADIOS2_CORE_LIBRARY}")
    endif()

    if(ADIOS2_INCLUDE_DIR AND ADIOS2_CORE_LIBRARY)
        set(_adios2_FOUND TRUE)
    endif()
endif()

# --- Standard find_package handling ---
if(_adios2_FOUND AND NOT ADIOS2_FOUND)
    # Only call this if we didn't already find it via config mode
    find_package_handle_standard_args(ADIOS2
        REQUIRED_VARS
            ADIOS2_LIBRARIES
            ADIOS2_INCLUDE_DIRS
        VERSION_VAR
            ADIOS2_VERSION
        FAIL_MESSAGE
            "Could not find ADIOS2. Set ADIOS2_DIR or ADIOS2_ROOT."
    )
elseif(NOT _adios2_FOUND)
    find_package_handle_standard_args(ADIOS2
        REQUIRED_VARS
            ADIOS2_LIBRARIES
            ADIOS2_INCLUDE_DIRS
        FAIL_MESSAGE
            "Could not find ADIOS2. Set ADIOS2_DIR or ADIOS2_ROOT."
    )
endif()

mark_as_advanced(
    ADIOS2_INCLUDE_DIR
    ADIOS2_CORE_LIBRARY
    ADIOS2_CORE_MPI_LIBRARY
    ADIOS2_FORTRAN_LIBRARY
    ADIOS2_FORTRAN_MPI_LIBRARY
)
