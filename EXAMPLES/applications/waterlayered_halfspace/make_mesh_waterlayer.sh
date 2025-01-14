#!/bin/bash

#
# Temporary instructions
#
# 1. set path to GEOCUBIT base directory (../../../CUBIT_GEOCUBIT/), for example:
#    export PYTHONPATH=$PYTHONPATH:/import/c/w/tape/3D/SPECFEM3D/CUBIT_GEOCUBIT
#    export PATH=$PATH:/import/c/w/tape/3D/SPECFEM3D/CUBIT_GEOCUBIT
#
#    check path:
#    which GEOCUBIT.py
#
# 2. run this script to generate mesh
#    ./make_mesh_waterlayer.sh
#

GEOCUBIT.py --build_volume --mesh --cfg=waterlayer_only.cfg
GEOCUBIT.py --collect --meshfiles=MESH_GEOCUBIT/mesh_vol_0.e --export2SPECFEM3D --SEMoutput=MESH
cp MESH/nummaterial_velocity_file.reference MESH/nummaterial_velocity_file
