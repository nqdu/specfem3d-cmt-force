# CompilerFlags.cmake
# --------------------
# Set up Fortran compiler flags based on the compiler ID and build type.
#
# This is the CMake equivalent of the flags.guess script used in the
# autotools-based SPECFEM3D build system.
#
# Variables set:
#   SPECFEM_Fortran_DEF_FLAGS   - Default flags (always applied)
#   SPECFEM_Fortran_OPT_FLAGS   - Optimization flags (Release builds)
#   SPECFEM_Fortran_DEBUG_FLAGS  - Debug flags (Debug builds)
#   SPECFEM_Fortran_OMP_FLAGS   - OpenMP flags
#   FLAGS_CHECK                  - Combined flags for the current build type
#
# Usage:
#   include(CompilerFlags)
#   # FLAGS_CHECK is automatically set based on CMAKE_BUILD_TYPE

include(CheckFortranCompilerFlag)

# ============================================================================
# Detect the real compiler behind potential wrappers (e.g., Cray ftn)
# ============================================================================
set(_specfem_FC_ID "${CMAKE_Fortran_COMPILER_ID}")

# On Cray systems, the ftn wrapper can hide the real compiler.
# Check the PE_ENV environment variable.
if(_specfem_FC_ID STREQUAL "Cray" OR
   CMAKE_Fortran_COMPILER MATCHES "ftn$")
    if(DEFINED ENV{PE_ENV})
        if("$ENV{PE_ENV}" STREQUAL "GNU")
            set(_specfem_FC_ID "GNU")
        elseif("$ENV{PE_ENV}" STREQUAL "INTEL")
            set(_specfem_FC_ID "Intel")
        elseif("$ENV{PE_ENV}" STREQUAL "PGI")
            set(_specfem_FC_ID "PGI")
        elseif("$ENV{PE_ENV}" STREQUAL "CRAY")
            set(_specfem_FC_ID "Cray")
        endif()
    endif()
endif()

# ============================================================================
# Default values
# ============================================================================
set(SPECFEM_Fortran_DEF_FLAGS "")
set(SPECFEM_Fortran_OPT_FLAGS "")
set(SPECFEM_Fortran_DEBUG_FLAGS "")
set(SPECFEM_Fortran_OMP_FLAGS "")

# ============================================================================
# Compiler-specific flags
# ============================================================================

if(_specfem_FC_ID STREQUAL "GNU")
    # -----------------------------------------------------------------------
    # GNU gfortran
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-std=f2008 -fimplicit-none -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-O3 -finline-functions"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS
        "-g -O0 -ggdb -fbacktrace -fbounds-check -frange-check -Werror"
    )
    set(SPECFEM_Fortran_OMP_FLAGS "-fopenmp")

elseif(_specfem_FC_ID STREQUAL "Intel" OR _specfem_FC_ID STREQUAL "IntelLLVM")
    # -----------------------------------------------------------------------
    # Intel ifort / ifx
    #
    # -assume buffered_io is important on parallel file systems (Lustre).
    # -xHost enables processor-specific optimizations.
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -std08 -diag-disable 6477 -implicitnone -gen-interfaces -warn all,noexternal"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-O3 -check nobounds"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS
        "-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv"
    )
    set(SPECFEM_Fortran_OMP_FLAGS "-qopenmp")

elseif(_specfem_FC_ID STREQUAL "PGI" OR _specfem_FC_ID STREQUAL "NVHPC")
    # -----------------------------------------------------------------------
    # Portland PGI / NVIDIA HPC SDK (pgfortran / nvfortran)
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-Mdclchk -Minform=warn -mcmodel=medium"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-Mnobounds -fast"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS
        "-Mbounds"
    )
    set(SPECFEM_Fortran_OMP_FLAGS "-mp")

elseif(_specfem_FC_ID STREQUAL "Cray")
    # -----------------------------------------------------------------------
    # Cray Fortran (crayftn)
    #
    # Aggressive optimization with Cray-specific flags.
    # OpenMP is enabled by default on Cray compilers.
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-M 1193 -M 1438"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-O3 -Onoaggress -Oipa0 -hfp2 -Ovector3 -Oscalar3 -Ocache2 -Ounroll2 -Ofusion2"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS
        "-eC -eD -ec -en -eI -ea -g -G0"
    )
    # OpenMP is enabled by default for Cray; empty flag
    set(SPECFEM_Fortran_OMP_FLAGS "")

elseif(_specfem_FC_ID STREQUAL "XL")
    # -----------------------------------------------------------------------
    # IBM XL Fortran (xlf)
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-qassert=contig -qhot -q64 -qtune=auto -qarch=auto -qcache=auto -qfree=f90 -qsuffix=f=f90 -qhalt=w -qlanglvl=2008std -qzerosize -g -qsuppress=1518-234 -qsuppress=1518-317 -qsuppress=1518-318 -qsuppress=1500-036 -qsuppress=1515-009"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-O4 -qstrict"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS
        "-g -O0 -C -qddim -qfullpath -qflttrap=overflow:zerodivide:invalid:enable -qfloat=nans -qinitauto=7FBFFFFF"
    )
    set(SPECFEM_Fortran_OMP_FLAGS "-qsmp=omp")

elseif(_specfem_FC_ID STREQUAL "PathScale")
    # -----------------------------------------------------------------------
    # PathScale pathf90
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS
        "-fno-math-errno -ffast-math -msse3 -march=auto -fno-second-underscore -align64"
    )
    set(SPECFEM_Fortran_OPT_FLAGS
        "-O3 -OPT:Ofast -LNO:fusion=2 -LNO:simd=2 -LNO:simd_verbose=ON"
    )
    set(SPECFEM_Fortran_DEBUG_FLAGS "-g2")
    set(SPECFEM_Fortran_OMP_FLAGS "-mp")

elseif(_specfem_FC_ID STREQUAL "LaheyFujitsu" OR
       CMAKE_Fortran_COMPILER MATCHES "lf95$")
    # -----------------------------------------------------------------------
    # Lahey/Fujitsu lf95
    # -----------------------------------------------------------------------
    set(SPECFEM_Fortran_DEF_FLAGS "--warn --wo --tpp --f95 --dal")
    set(SPECFEM_Fortran_OPT_FLAGS "-O")
    set(SPECFEM_Fortran_DEBUG_FLAGS "--chk")
    set(SPECFEM_Fortran_OMP_FLAGS "")

else()
    # -----------------------------------------------------------------------
    # Unknown compiler - use conservative defaults
    # -----------------------------------------------------------------------
    message(STATUS "CompilerFlags: Unrecognized Fortran compiler '${_specfem_FC_ID}'. Using conservative defaults.")
    set(SPECFEM_Fortran_DEF_FLAGS "")
    set(SPECFEM_Fortran_OPT_FLAGS "-O2")
    set(SPECFEM_Fortran_DEBUG_FLAGS "-g -O0")
    set(SPECFEM_Fortran_OMP_FLAGS "")
endif()

# ============================================================================
# Set FLAGS_CHECK based on CMAKE_BUILD_TYPE
# ============================================================================
# The FLAGS_CHECK variable mirrors the autotools build system convention.
# It combines the default flags with either optimization or debug flags.
#
# Build types recognized:
#   Debug, RelWithDebInfo -> use debug flags
#   Release, MinSizeRel, "" (default) -> use optimization flags

if(NOT DEFINED FLAGS_CHECK OR FLAGS_CHECK STREQUAL "")
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type_upper)

    if(_build_type_upper STREQUAL "DEBUG" OR _build_type_upper STREQUAL "RELWITHDEBINFO")
        set(FLAGS_CHECK "${SPECFEM_Fortran_DEF_FLAGS} ${SPECFEM_Fortran_DEBUG_FLAGS}"
            CACHE STRING "Fortran compiler flags for SPECFEM3D" FORCE)
        message(STATUS "CompilerFlags: Using DEBUG flags for ${_specfem_FC_ID}")
    else()
        set(FLAGS_CHECK "${SPECFEM_Fortran_DEF_FLAGS} ${SPECFEM_Fortran_OPT_FLAGS}"
            CACHE STRING "Fortran compiler flags for SPECFEM3D" FORCE)
        message(STATUS "CompilerFlags: Using RELEASE/OPT flags for ${_specfem_FC_ID}")
    endif()
endif()

# ============================================================================
# OpenMP flags
# ============================================================================
if(NOT DEFINED OMP_FCFLAGS OR OMP_FCFLAGS STREQUAL "")
    set(OMP_FCFLAGS "${SPECFEM_Fortran_OMP_FLAGS}"
        CACHE STRING "OpenMP Fortran flags for SPECFEM3D")
endif()

# ============================================================================
# Report
# ============================================================================
message(STATUS "CompilerFlags: Fortran compiler ID = ${_specfem_FC_ID}")
message(STATUS "CompilerFlags: FLAGS_CHECK = ${FLAGS_CHECK}")
if(SPECFEM_Fortran_OMP_FLAGS)
    message(STATUS "CompilerFlags: OMP_FCFLAGS = ${OMP_FCFLAGS}")
endif()
