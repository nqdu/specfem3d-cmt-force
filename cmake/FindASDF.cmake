# FindASDF.cmake
# ---------------
# Find the ASDF (Adaptable Seismic Data Format) library.
#
# ASDF provides HDF5-based seismic data storage, commonly used with
# SPECFEM3D for waveform output.
#
# Hints:
#   ASDF_DIR       - Root directory of ASDF installation (env or CMake var)
#   ASDF_INCLUDEDIR - Directory containing ASDF headers/modules
#   ASDF_LIBDIR     - Directory containing ASDF libraries
#
# Result variables:
#   ASDF_FOUND         - True if ASDF was found
#   ASDF_INCLUDE_DIRS  - Include directories for ASDF (headers and Fortran modules)
#   ASDF_LIBRARIES     - Libraries to link against

include(FindPackageHandleStandardArgs)

# Build search path list
set(_asdf_SEARCH_DIRS "")
if(ASDF_DIR)
    list(APPEND _asdf_SEARCH_DIRS "${ASDF_DIR}")
endif()
if(DEFINED ENV{ASDF_DIR})
    list(APPEND _asdf_SEARCH_DIRS "$ENV{ASDF_DIR}")
endif()

# --- Find ASDF header or Fortran module ---
# ASDF can be a C library (asdf.h) or a Fortran library (with .mod files)
find_path(ASDF_INCLUDE_DIR
    NAMES asdf.h ASDF.h asdf_data.mod
    HINTS
        ${ASDF_INCLUDEDIR}
        $ENV{ASDF_INCLUDEDIR}
    PATHS
        ${_asdf_SEARCH_DIRS}
    PATH_SUFFIXES
        include
    DOC "Directory containing ASDF headers or Fortran module files"
)

# --- Find libasdf ---
find_library(ASDF_LIBRARY
    NAMES asdf
    HINTS
        ${ASDF_LIBDIR}
        $ENV{ASDF_LIBDIR}
    PATHS
        ${_asdf_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The ASDF library"
)

# --- Find libasdf_fortran (some installations split Fortran bindings) ---
find_library(ASDF_FORTRAN_LIBRARY
    NAMES asdf_fortran asdf-fortran
    HINTS
        ${ASDF_LIBDIR}
        $ENV{ASDF_LIBDIR}
    PATHS
        ${_asdf_SEARCH_DIRS}
    PATH_SUFFIXES
        lib
        lib64
    DOC "The ASDF Fortran bindings library"
)

# --- Collect results ---
set(ASDF_INCLUDE_DIRS "")
if(ASDF_INCLUDE_DIR)
    set(ASDF_INCLUDE_DIRS "${ASDF_INCLUDE_DIR}")
endif()

set(ASDF_LIBRARIES "")
if(ASDF_FORTRAN_LIBRARY)
    list(APPEND ASDF_LIBRARIES "${ASDF_FORTRAN_LIBRARY}")
endif()
if(ASDF_LIBRARY)
    list(APPEND ASDF_LIBRARIES "${ASDF_LIBRARY}")
endif()

# --- Standard find_package handling ---
find_package_handle_standard_args(ASDF
    REQUIRED_VARS
        ASDF_LIBRARY
    FAIL_MESSAGE
        "Could not find ASDF library. Set ASDF_DIR or ASDF_INCLUDEDIR/ASDF_LIBDIR."
)

mark_as_advanced(
    ASDF_INCLUDE_DIR
    ASDF_LIBRARY
    ASDF_FORTRAN_LIBRARY
)
