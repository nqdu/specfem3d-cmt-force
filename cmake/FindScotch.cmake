# FindScotch.cmake
# -----------------
# Find the PT-Scotch and Scotch graph partitioning libraries.
#
# This module looks for both the parallel (PT-Scotch) and sequential (Scotch)
# libraries and their associated headers.
#
# Hints:
#   SCOTCH_DIR       - Root directory of Scotch installation (env or CMake var)
#   SCOTCH_INCLUDEDIR - Directory containing Scotch headers
#   SCOTCH_LIBDIR     - Directory containing Scotch libraries
#
# Result variables:
#   SCOTCH_FOUND            - True if Scotch was found
#   SCOTCH_INCLUDE_DIRS     - Include directories for Scotch
#   SCOTCH_LIBRARIES        - Libraries to link against
#   SCOTCH_HAS_PTSCOTCH     - True if PT-Scotch (parallel) was found
#
# Components:
#   scotch, scotcherr, ptscotch, ptscotcherr

include(FindPackageHandleStandardArgs)

# Build a list of search paths from hints
set(_scotch_SEARCH_DIRS "")

if(SCOTCH_DIR)
    list(APPEND _scotch_SEARCH_DIRS "${SCOTCH_DIR}")
endif()

if(DEFINED ENV{SCOTCH_DIR})
    list(APPEND _scotch_SEARCH_DIRS "$ENV{SCOTCH_DIR}")
endif()

# --- Find scotch.h ---
find_path(SCOTCH_INCLUDE_DIR
    NAMES scotch.h
    HINTS
        ${SCOTCH_INCLUDEDIR}
        $ENV{SCOTCH_INCLUDEDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        include
        include/scotch
    DOC "Directory containing scotch.h"
)

# --- Find ptscotch.h ---
find_path(PTSCOTCH_INCLUDE_DIR
    NAMES ptscotch.h
    HINTS
        ${SCOTCH_INCLUDEDIR}
        $ENV{SCOTCH_INCLUDEDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        include
        include/scotch
    DOC "Directory containing ptscotch.h"
)

# --- Find libscotch ---
find_library(SCOTCH_LIBRARY
    NAMES scotch
    HINTS
        ${SCOTCH_LIBDIR}
        $ENV{SCOTCH_LIBDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The Scotch library"
)

# --- Find libscotcherr ---
find_library(SCOTCHERR_LIBRARY
    NAMES scotcherr
    HINTS
        ${SCOTCH_LIBDIR}
        $ENV{SCOTCH_LIBDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The Scotch error-handling library"
)

# --- Find libptscotch ---
find_library(PTSCOTCH_LIBRARY
    NAMES ptscotch
    HINTS
        ${SCOTCH_LIBDIR}
        $ENV{SCOTCH_LIBDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The PT-Scotch (parallel) library"
)

# --- Find libptscotcherr ---
find_library(PTSCOTCHERR_LIBRARY
    NAMES ptscotcherr
    HINTS
        ${SCOTCH_LIBDIR}
        $ENV{SCOTCH_LIBDIR}
    PATHS
        ${_scotch_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The PT-Scotch error-handling library"
)

# --- Determine if PT-Scotch is available ---
set(SCOTCH_HAS_PTSCOTCH FALSE)
if(PTSCOTCH_INCLUDE_DIR AND PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY)
    set(SCOTCH_HAS_PTSCOTCH TRUE)
endif()

# --- Collect include directories ---
set(SCOTCH_INCLUDE_DIRS "")
if(SCOTCH_INCLUDE_DIR)
    list(APPEND SCOTCH_INCLUDE_DIRS "${SCOTCH_INCLUDE_DIR}")
endif()
if(PTSCOTCH_INCLUDE_DIR AND NOT "${PTSCOTCH_INCLUDE_DIR}" STREQUAL "${SCOTCH_INCLUDE_DIR}")
    list(APPEND SCOTCH_INCLUDE_DIRS "${PTSCOTCH_INCLUDE_DIR}")
endif()

# --- Collect libraries ---
# Order matters: ptscotch before scotch for proper linking
set(SCOTCH_LIBRARIES "")
if(SCOTCH_HAS_PTSCOTCH)
    list(APPEND SCOTCH_LIBRARIES "${PTSCOTCH_LIBRARY}" "${PTSCOTCHERR_LIBRARY}")
endif()
if(SCOTCH_LIBRARY)
    list(APPEND SCOTCH_LIBRARIES "${SCOTCH_LIBRARY}")
endif()
if(SCOTCHERR_LIBRARY)
    list(APPEND SCOTCH_LIBRARIES "${SCOTCHERR_LIBRARY}")
endif()

# --- Standard find_package handling ---
find_package_handle_standard_args(Scotch
    REQUIRED_VARS
        SCOTCH_LIBRARY
        SCOTCHERR_LIBRARY
        SCOTCH_INCLUDE_DIR
    FAIL_MESSAGE
        "Could not find Scotch. Set SCOTCH_DIR or SCOTCH_INCLUDEDIR/SCOTCH_LIBDIR."
)

mark_as_advanced(
    SCOTCH_INCLUDE_DIR
    PTSCOTCH_INCLUDE_DIR
    SCOTCH_LIBRARY
    SCOTCHERR_LIBRARY
    PTSCOTCH_LIBRARY
    PTSCOTCHERR_LIBRARY
)
