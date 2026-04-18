# FindMETIS.cmake
# -----------------
# Find the METIS graph partitioning library.
#
# METIS is used for mesh partitioning/decomposition in SPECFEM3D.
#
# Hints:
#   METIS_DIR        - Root directory of METIS installation (env or CMake var)
#   METIS_INCLUDEDIR - Directory containing metis.h
#   METIS_LIBDIR     - Directory containing METIS libraries
#
# Result variables:
#   METIS_FOUND         - True if METIS was found
#   METIS_INCLUDE_DIRS  - Include directories for METIS
#   METIS_LIBRARIES     - Libraries to link against
#   METIS_VERSION       - METIS version string (if determinable from header)

include(FindPackageHandleStandardArgs)

# Build search path list
set(_metis_SEARCH_DIRS "")
if(METIS_DIR)
    list(APPEND _metis_SEARCH_DIRS "${METIS_DIR}")
endif()
if(DEFINED ENV{METIS_DIR})
    list(APPEND _metis_SEARCH_DIRS "$ENV{METIS_DIR}")
endif()

# --- Find metis.h ---
find_path(METIS_INCLUDE_DIR
    NAMES metis.h
    HINTS
        ${METIS_INCLUDEDIR}
        $ENV{METIS_INCLUDEDIR}
    PATHS
        ${_metis_SEARCH_DIRS}
    PATH_SUFFIXES
        include
    DOC "Directory containing metis.h"
)

# --- Find libmetis ---
find_library(METIS_LIBRARY
    NAMES metis
    HINTS
        ${METIS_LIBDIR}
        $ENV{METIS_LIBDIR}
    PATHS
        ${_metis_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The METIS library"
)

# --- Try to determine version from metis.h ---
if(METIS_INCLUDE_DIR AND EXISTS "${METIS_INCLUDE_DIR}/metis.h")
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_version_major
        REGEX "^#define[ \t]+METIS_VER_MAJOR[ \t]+[0-9]+")
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_version_minor
        REGEX "^#define[ \t]+METIS_VER_MINOR[ \t]+[0-9]+")
    file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_version_subminor
        REGEX "^#define[ \t]+METIS_VER_SUBMINOR[ \t]+[0-9]+")

    if(_metis_version_major)
        string(REGEX REPLACE "^#define[ \t]+METIS_VER_MAJOR[ \t]+([0-9]+)" "\\1"
            METIS_VERSION_MAJOR "${_metis_version_major}")
        string(REGEX REPLACE "^#define[ \t]+METIS_VER_MINOR[ \t]+([0-9]+)" "\\1"
            METIS_VERSION_MINOR "${_metis_version_minor}")
        string(REGEX REPLACE "^#define[ \t]+METIS_VER_SUBMINOR[ \t]+([0-9]+)" "\\1"
            METIS_VERSION_SUBMINOR "${_metis_version_subminor}")
        set(METIS_VERSION "${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_SUBMINOR}")
    endif()
endif()

# --- Collect results ---
set(METIS_INCLUDE_DIRS "")
if(METIS_INCLUDE_DIR)
    set(METIS_INCLUDE_DIRS "${METIS_INCLUDE_DIR}")
endif()

set(METIS_LIBRARIES "")
if(METIS_LIBRARY)
    set(METIS_LIBRARIES "${METIS_LIBRARY}")
endif()

# --- Standard find_package handling ---
find_package_handle_standard_args(METIS
    REQUIRED_VARS
        METIS_LIBRARY
        METIS_INCLUDE_DIR
    VERSION_VAR
        METIS_VERSION
    FAIL_MESSAGE
        "Could not find METIS. Set METIS_DIR or METIS_INCLUDEDIR/METIS_LIBDIR."
)

mark_as_advanced(
    METIS_INCLUDE_DIR
    METIS_LIBRARY
)
