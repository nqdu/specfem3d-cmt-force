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

# static_gpu.cu: only built (against src/gpu's CUDA/HIP headers) when GPU support is on,
# otherwise static_module.f90 links against the no-op stubs in
# src/gpu/specfem3D_gpu_cuda_method_stubs.c (part of specfem3D_SPECFEM_OBJECTS already)
ifeq ($(HAS_GPU),yes)
  ifeq ($(CUDA),yes)
    # NVCC_FLAGS uses -dc (relocatable device code), so the raw compiled object still needs
    # its own "nvcc -dlink" pass before the host linker can use it -- mirrors how
    # cuda_specfem3D_DEVICE_OBJ is built from gpu_specfem3D_OBJECTS in src/gpu/rules.mk
    static3D_GPU_OBJECT = $O/static_gpu.static_cuda.o $O/static_gpu_dlink.static_cuda.o
  endif
  ifeq ($(HIP),yes)
    static3D_GPU_OBJECT = $O/static_gpu.static_hip.o
  endif
else
static3D_GPU_OBJECT =
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
	$O/static_init.static_module.o \
	$O/static_impl.static_module.o \
	$O/xstatic3D.static.o \
	$(static3D_PETSC_OBJECT) \
	$(static3D_GPU_OBJECT) \
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

## submodules (static_impl.f90/static_init.f90) implement the module procedures declared in
## static_module.f90, so they must be compiled after it produces its .mod/.smod
$O/static_init.static_module.o: $O/static_module.static_module.o
$O/static_impl.static_module.o: $O/static_module.static_module.o

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

## static_gpu.cu: uses src/gpu's mesh_constants_gpu.h / CUDA / HIP headers and NVCC_FLAGS,
## SELECTOR_CFLAG etc. defined in src/gpu/rules.mk; mesh_constants_cuda.h pulls in
## kernels/kernel_proto.cu.h, so the kernels/ dir needs to be on the include path too
GPU_S := ${S_TOP}/src/gpu
GPU_KERNEL_DIR := $(GPU_S)/kernels

$O/static_gpu.static_cuda.o: $S/static_gpu.cu ${SETUP}/config.h $(GPU_S)/mesh_constants_gpu.h $(GPU_S)/mesh_constants_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS) -I${SETUP} -I$(GPU_S) -I$(GPU_KERNEL_DIR) $(SELECTOR_CFLAG)

$O/static_gpu_dlink.static_cuda.o: $O/static_gpu.static_cuda.o
	${NVCCLINK} -o $@ $^

$O/static_gpu.static_hip.o: $S/static_gpu.cu ${SETUP}/config.h $(GPU_S)/mesh_constants_gpu.h $(GPU_S)/mesh_constants_hip.h
	${HIPCC} ${HIP_CFLAG_ENDING} -c $< -o $@ $(HIPCC_CFLAGS) -I${SETUP} -I$(GPU_S) -I$(GPU_KERNEL_DIR) $(SELECTOR_CFLAG)
