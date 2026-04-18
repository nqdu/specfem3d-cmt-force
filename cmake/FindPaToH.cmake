# FindPaToH.cmake
# -----------------
# Find the PaToH (Partitioning Tool for Hypergraphs) library.
#
# PaToH is a hypergraph partitioning tool that can be used for mesh
# decomposition in SPECFEM3D.
#
# Hints:
#   PATOH_DIR        - Root directory of PaToH installation (env or CMake var)
#   PATOH_INCLUDEDIR - Directory containing patoh.h
#   PATOH_LIBDIR     - Directory containing PaToH libraries
#
# Result variables:
#   PATOH_FOUND         - True if PaToH was found
#   PATOH_INCLUDE_DIRS  - Include directories for PaToH
#   PATOH_LIBRARIES     - Libraries to link against

include(FindPackageHandleStandardArgs)

# Build search path list
set(_patoh_SEARCH_DIRS "")
if(PATOH_DIR)
    list(APPEND _patoh_SEARCH_DIRS "${PATOH_DIR}")
endif()
if(DEFINED ENV{PATOH_DIR})
    list(APPEND _patoh_SEARCH_DIRS "$ENV{PATOH_DIR}")
endif()

# --- Find patoh.h ---
find_path(PATOH_INCLUDE_DIR
    NAMES patoh.h
    HINTS
        ${PATOH_INCLUDEDIR}
        $ENV{PATOH_INCLUDEDIR}
    PATHS
        ${_patoh_SEARCH_DIRS}
    PATH_SUFFIXES
        include
    DOC "Directory containing patoh.h"
)

# --- Find libpatoh ---
find_library(PATOH_LIBRARY
    NAMES patoh
    HINTS
        ${PATOH_LIBDIR}
        $ENV{PATOH_LIBDIR}
    PATHS
        ${_patoh_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The PaToH library"
)

# --- Collect results ---
set(PATOH_INCLUDE_DIRS "")
if(PATOH_INCLUDE_DIR)
    set(PATOH_INCLUDE_DIRS "${PATOH_INCLUDE_DIR}")
endif()

set(PATOH_LIBRARIES "")
if(PATOH_LIBRARY)
    set(PATOH_LIBRARIES "${PATOH_LIBRARY}")
endif()

# --- Standard find_package handling ---
find_package_handle_standard_args(PaToH
    REQUIRED_VARS
        PATOH_LIBRARY
        PATOH_INCLUDE_DIR
    FAIL_MESSAGE
        "Could not find PaToH. Set PATOH_DIR or PATOH_INCLUDEDIR/PATOH_LIBDIR."
)

# PaToH often needs libm -- must be after find_package_handle_standard_args sets PaToH_FOUND
if(PaToH_FOUND)
    find_library(_patoh_math_lib NAMES m)
    if(_patoh_math_lib)
        list(APPEND PATOH_LIBRARIES "${_patoh_math_lib}")
    endif()
endif()

mark_as_advanced(
    PATOH_INCLUDE_DIR
    PATOH_LIBRARY
)
