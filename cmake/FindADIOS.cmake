# FindADIOS.cmake
# ----------------
# Find the ADIOS 1.x I/O library.
#
# ADIOS 1.x provides an `adios_config` tool that returns compiler and linker
# flags. This module tries that first, then falls back to manual header/library
# search.
#
# Hints:
#   ADIOS_DIR       - Root directory of ADIOS installation (env or CMake var)
#   ADIOS_INCLUDEDIR - Directory containing ADIOS headers
#   ADIOS_LIBDIR     - Directory containing ADIOS libraries
#
# Result variables:
#   ADIOS_FOUND         - True if ADIOS was found
#   ADIOS_INCLUDE_DIRS  - Include directories for ADIOS
#   ADIOS_LIBRARIES     - Libraries to link against
#   ADIOS_VERSION       - ADIOS version string (if available)

include(FindPackageHandleStandardArgs)

set(_adios_FOUND_VIA_CONFIG FALSE)

# Build search path list
set(_adios_SEARCH_DIRS "")
if(ADIOS_DIR)
    list(APPEND _adios_SEARCH_DIRS "${ADIOS_DIR}")
endif()
if(DEFINED ENV{ADIOS_DIR})
    list(APPEND _adios_SEARCH_DIRS "$ENV{ADIOS_DIR}")
endif()

# --- Strategy 1: Use adios_config tool ---
# Look for adios_config in PATH or under ADIOS_DIR/bin
find_program(ADIOS_CONFIG_EXECUTABLE
    NAMES adios_config
    HINTS
        ${_adios_SEARCH_DIRS}
    PATH_SUFFIXES
        bin
    DOC "ADIOS 1.x configuration tool"
)

if(ADIOS_CONFIG_EXECUTABLE)
    # Get version
    execute_process(
        COMMAND ${ADIOS_CONFIG_EXECUTABLE} -v
        OUTPUT_VARIABLE _adios_version_output
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _adios_config_version_result
    )
    if(_adios_config_version_result EQUAL 0)
        set(ADIOS_VERSION "${_adios_version_output}")
    endif()

    # Get C compiler flags (includes)
    execute_process(
        COMMAND ${ADIOS_CONFIG_EXECUTABLE} -c
        OUTPUT_VARIABLE _adios_cflags
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _adios_config_c_result
    )

    # Get Fortran flags
    execute_process(
        COMMAND ${ADIOS_CONFIG_EXECUTABLE} -f -c
        OUTPUT_VARIABLE _adios_fflags
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _adios_config_fc_result
    )

    # Get linker flags
    execute_process(
        COMMAND ${ADIOS_CONFIG_EXECUTABLE} -l
        OUTPUT_VARIABLE _adios_ldflags
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _adios_config_l_result
    )

    # Get Fortran linker flags
    execute_process(
        COMMAND ${ADIOS_CONFIG_EXECUTABLE} -f -l
        OUTPUT_VARIABLE _adios_fldflags
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _adios_config_fl_result
    )

    # Parse include directories from -I flags
    set(ADIOS_INCLUDE_DIRS "")
    if(_adios_config_c_result EQUAL 0)
        string(REGEX MATCHALL "-I[^ ]+" _adios_inc_flags "${_adios_cflags}")
        foreach(_flag ${_adios_inc_flags})
            string(REGEX REPLACE "^-I" "" _dir "${_flag}")
            list(APPEND ADIOS_INCLUDE_DIRS "${_dir}")
        endforeach()
    endif()
    # Also parse Fortran include dirs
    if(_adios_config_fc_result EQUAL 0)
        string(REGEX MATCHALL "-I[^ ]+" _adios_finc_flags "${_adios_fflags}")
        foreach(_flag ${_adios_finc_flags})
            string(REGEX REPLACE "^-I" "" _dir "${_flag}")
            list(APPEND ADIOS_INCLUDE_DIRS "${_dir}")
        endforeach()
    endif()
    if(ADIOS_INCLUDE_DIRS)
        list(REMOVE_DUPLICATES ADIOS_INCLUDE_DIRS)
    endif()

    # Parse libraries from linker flags
    set(ADIOS_LIBRARIES "")
    set(_adios_all_ldflags "")
    if(_adios_config_fl_result EQUAL 0)
        set(_adios_all_ldflags "${_adios_fldflags}")
    elseif(_adios_config_l_result EQUAL 0)
        set(_adios_all_ldflags "${_adios_ldflags}")
    endif()

    if(_adios_all_ldflags)
        # Extract -L paths
        set(_adios_lib_dirs "")
        string(REGEX MATCHALL "-L[^ ]+" _adios_libdir_flags "${_adios_all_ldflags}")
        foreach(_flag ${_adios_libdir_flags})
            string(REGEX REPLACE "^-L" "" _dir "${_flag}")
            list(APPEND _adios_lib_dirs "${_dir}")
        endforeach()

        # Extract -l libraries
        string(REGEX MATCHALL "-l[^ ]+" _adios_lib_flags "${_adios_all_ldflags}")
        foreach(_flag ${_adios_lib_flags})
            string(REGEX REPLACE "^-l" "" _lib "${_flag}")
            find_library(_adios_lib_${_lib}
                NAMES ${_lib}
                HINTS ${_adios_lib_dirs}
            )
            if(_adios_lib_${_lib})
                list(APPEND ADIOS_LIBRARIES "${_adios_lib_${_lib}}")
            else()
                # Fall back to just the -l flag
                list(APPEND ADIOS_LIBRARIES "${_flag}")
            endif()
        endforeach()

        # Also include any raw .a or .so files on the link line
        string(REGEX MATCHALL "[^ ]+\\.(a|so|dylib)" _adios_raw_libs "${_adios_all_ldflags}")
        foreach(_lib ${_adios_raw_libs})
            list(APPEND ADIOS_LIBRARIES "${_lib}")
        endforeach()
    endif()

    if(ADIOS_INCLUDE_DIRS AND ADIOS_LIBRARIES)
        set(_adios_FOUND_VIA_CONFIG TRUE)
    endif()
endif()

# --- Strategy 2: Manual search (fallback) ---
if(NOT _adios_FOUND_VIA_CONFIG)
    find_path(ADIOS_INCLUDE_DIR
        NAMES adios.h
        HINTS
            ${ADIOS_INCLUDEDIR}
            $ENV{ADIOS_INCLUDEDIR}
        PATHS
            ${_adios_SEARCH_DIRS}
        PATH_SUFFIXES
            include
        DOC "Directory containing adios.h"
    )

    find_library(ADIOS_LIBRARY
        NAMES adios
        HINTS
            ${ADIOS_LIBDIR}
            $ENV{ADIOS_LIBDIR}
        PATHS
            ${_adios_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS library"
    )

    # ADIOS 1.x also needs adiosf for Fortran bindings
    find_library(ADIOSF_LIBRARY
        NAMES adiosf
        HINTS
            ${ADIOS_LIBDIR}
            $ENV{ADIOS_LIBDIR}
        PATHS
            ${_adios_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS Fortran library"
    )

    # Read-only transport (often needed)
    find_library(ADIOSREAD_LIBRARY
        NAMES adiosread
        HINTS
            ${ADIOS_LIBDIR}
            $ENV{ADIOS_LIBDIR}
        PATHS
            ${_adios_SEARCH_DIRS}
        PATH_SUFFIXES
            lib
            lib64
        DOC "The ADIOS read library"
    )

    set(ADIOS_INCLUDE_DIRS "")
    if(ADIOS_INCLUDE_DIR)
        set(ADIOS_INCLUDE_DIRS "${ADIOS_INCLUDE_DIR}")
    endif()

    set(ADIOS_LIBRARIES "")
    if(ADIOSF_LIBRARY)
        list(APPEND ADIOS_LIBRARIES "${ADIOSF_LIBRARY}")
    endif()
    if(ADIOS_LIBRARY)
        list(APPEND ADIOS_LIBRARIES "${ADIOS_LIBRARY}")
    endif()
    if(ADIOSREAD_LIBRARY)
        list(APPEND ADIOS_LIBRARIES "${ADIOSREAD_LIBRARY}")
    endif()
endif()

# --- Standard find_package handling ---
if(_adios_FOUND_VIA_CONFIG)
    find_package_handle_standard_args(ADIOS
        REQUIRED_VARS
            ADIOS_LIBRARIES
            ADIOS_INCLUDE_DIRS
        VERSION_VAR
            ADIOS_VERSION
        FAIL_MESSAGE
            "Could not find ADIOS 1.x. Set ADIOS_DIR or ensure adios_config is in PATH."
    )
else()
    find_package_handle_standard_args(ADIOS
        REQUIRED_VARS
            ADIOS_LIBRARY
            ADIOS_INCLUDE_DIR
        FAIL_MESSAGE
            "Could not find ADIOS 1.x. Set ADIOS_DIR or ADIOS_INCLUDEDIR/ADIOS_LIBDIR."
    )
endif()

mark_as_advanced(
    ADIOS_CONFIG_EXECUTABLE
    ADIOS_INCLUDE_DIR
    ADIOS_LIBRARY
    ADIOSF_LIBRARY
    ADIOSREAD_LIBRARY
)
