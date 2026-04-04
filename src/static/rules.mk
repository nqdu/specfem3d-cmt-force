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

#######################################

####
#### targets
####

static3D_TARGETS = \
	$E/xstatic3D \
	$(EMPTY_MACRO)

static3D_OBJECTS = \
	$O/static_module.static_module.o \
	$O/xstatic3D.static.o \
	$(EMPTY_MACRO)

# specfem3D objects without the main program entry point
static3D_SPECFEM_OBJECTS = $(filter-out $O/specfem3D.spec.o, $(specfem3D_OBJECTS))

static3D_MODULES = \
	$(FC_MODDIR)/static_module.$(FC_MODEXT) \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####

static3D: xstatic3D
xstatic3D: $E/xstatic3D

$E/xstatic3D: $(static3D_OBJECTS) $(static3D_SPECFEM_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo "building xstatic3D"
	@echo ""
	${FCLINK} -o $@ $(static3D_OBJECTS) $(static3D_SPECFEM_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(MPILIBS) $(VTKLIBS) $(SPECFEM_LINK_FLAGS)
	@echo ""

#######################################

####
#### rule to build each .o file below
####

## module file: depends on specfem3D_par and pml_par modules
$O/%.static_module.o: $S/%.f90 $O/specfem3D_par.spec_module.o $O/pml_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## main program: depends on static_module and specfem3D_par module
$O/%.static.o: $S/%.f90 $O/static_module.static_module.o $O/specfem3D_par.spec_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
