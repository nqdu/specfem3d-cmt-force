#=====================================================================
#
#                         S p e c f e m 3 D
#                         -----------------
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

## compilation directories
S := ${S_TOP}/src/static
$(static3D_OBJECTS): S = ${S_TOP}/src/static

ifneq ($(filter 1 yes YES true TRUE on ON,$(WITH_PETSC)),)
static3D_PETSC_OBJECT = $O/petsc_routines.static_c.o
static3D_PETSC_INCLUDE = $(PETSC_CPPFLAGS)
static3D_EXTRA_LIBS = $(PETSCLIBS)
static3D_CPPFLAGS = $(FC_DEFINE)WITH_PETSC $(PETSC_CPPFLAGS)
else
static3D_PETSC_OBJECT = $O/petsc_routines_stubs.static_c.o
static3D_PETSC_INCLUDE =
static3D_EXTRA_LIBS =
static3D_CPPFLAGS =
endif

#######################################

####
#### targets
####

static3D_TARGETS = \
	$E/xstatic3D \
	$(EMPTY_MACRO)

static3D_OBJECTS = \
	$O/petsc_interfaces.static_module.o \
	$O/static_module.static_module.o \
	$O/xstatic3D.static.o \
	$(static3D_PETSC_OBJECT) \
	$(EMPTY_MACRO)

# specfem3D objects without the main program entry point
static3D_SPECFEM_OBJECTS = $(filter-out $O/specfem3D.spec.o, $(specfem3D_OBJECTS))

static3D_MODULES = \
	$(FC_MODDIR)/petsc_interfaces.$(FC_MODEXT) \
	$(FC_MODDIR)/static_module.$(FC_MODEXT) \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####

static3D_LIBS = $(MPILIBS) $(VTKLIBS) $(static3D_EXTRA_LIBS)

ifeq ($(HAS_GPU),yes)
static3D_LIBS += $(GPU_LINK)
INFO_STATIC3D = "building xstatic3D $(BUILD_VERSION_TXT)"
else
INFO_STATIC3D = "building xstatic3D"
endif

static3D: xstatic3D
xstatic3D: $E/xstatic3D

$E/xstatic3D: $(static3D_OBJECTS) $(static3D_SPECFEM_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_STATIC3D)
	@echo ""
	${FCLINK} -o $@ $(static3D_OBJECTS) $(static3D_SPECFEM_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(static3D_LIBS) $(SPECFEM_LINK_FLAGS)
	@echo ""

#######################################

####
#### rule to build each .o file below
####

$O/static_module.static_module.o: $O/petsc_interfaces.static_module.o

$O/petsc_interfaces.static_module.o: $S/petsc_interfaces.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(static3D_CPPFLAGS) -c -o $@ $<

## module file: depends on specfem3D_par and pml_par modules
$O/%.static_module.o: $S/%.f90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(static3D_CPPFLAGS) -c -o $@ $<

## main program: depends on static_module and specfem3D_par module
$O/%.static.o: $S/%.f90 $O/static_module.static_module.o $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(static3D_CPPFLAGS) -c -o $@ $<

$O/%.static.o: $S/%.F90 $O/static_module.static_module.o $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(static3D_CPPFLAGS) -c -o $@ $<

$O/petsc_routines.static_c.o: $S/petsc_routines.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(static3D_PETSC_INCLUDE) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

$O/petsc_routines_stubs.static_c.o: $S/petsc_routines_stubs.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<
