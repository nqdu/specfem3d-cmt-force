# SPECFEM3D Cartesian

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10823181.svg)](https://doi.org/10.5281/zenodo.10823181)

SPECFEM3D Cartesian simulates acoustic (fluid), elastic (solid), coupled acoustic/elastic, poroelastic or seismic wave propagation in any type of conforming mesh of hexahedra (structured or not.)

It can, for instance, model seismic waves propagating in sedimentary basins or any other regional geological model following earthquakes. It can also be used for non-destructive testing or for ocean acoustics


SPECFEM3D Cartesian was founded by Dimitri Komatitsch and Jeroen Tromp, and is now being developed by a large, collaborative, and inclusive community. A complete list of authors can be found at
https://specfem3d.readthedocs.io/en/latest/authors/


## Installation

Instructions on how to install and use SPECFEM3D Cartesian are
available in the

- PDF manual located in directory: [doc/USER_MANUAL](doc/USER_MANUAL)

- HTML manual (latest version): [specfem3d.readthedocs.io](http://specfem3d.readthedocs.io/)


For a quick test, run the default example with these commands:
```
./configure FC=gfortran CC=gcc
make all
./run_this_example.sh
```
and check the output files in `./OUTPUT_FILES/`

>__NOTE__: Do not modify the 'configure' script directly. Please modify the
    'configure.ac' file instead, and generate a new 'configure' script with
    the command: `autoreconf -i`


## Development

[![Actions Status](https://github.com/SPECFEM/specfem3d/workflows/CI/badge.svg)](https://github.com/SPECFEM/specfem3d/actions)
[![Travis Status](https://app.travis-ci.com/SPECFEM/specfem3d.svg?branch=devel)](https://app.travis-ci.com/SPECFEM/specfem3d)
[![Azure Status](https://dev.azure.com/danielpeter22/SPECFEM3D/_apis/build/status/geodynamics.specfem3d?branchName=devel)](https://dev.azure.com/danielpeter22/SPECFEM3D/_build/latest?definitionId=5&branchName=devel)
[![codecov](https://codecov.io/gh/SPECFEM/specfem3d/branch/devel/graph/badge.svg)](https://codecov.io/gh/SPECFEM/specfem3d)
[![Documentation Status](https://readthedocs.org/projects/specfem3d/badge/?version=latest)](https://specfem3d.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

* Actions tests: [github actions specfem3d](https://github.com/SPECFEM/specfem3d/actions)

* Travis tests: [travis-ci specfem3d](https://travis-ci.com/SPECFEM/specfem3d/builds)


Development is hosted on GitHub in the
[SPECFEM/specfem3d repository](https://github.com/SPECFEM/specfem3d).

SPECFEM3D Cartesian is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive and
participatory community so we are happy that you are interested in
participating! We have collected a set of guidelines and advice on how to get
involved in the community and keep them in the specfem3d github wiki:
[specfem3d wiki](https://github.com/SPECFEM/specfem3d/wiki)


## Our contributors :sparkles:

[![SPECFEM3D Cartesian contributors](https://contrib.rocks/image?repo=SPECFEM/specfem3d)](https://github.com/SPECFEM/specfem3d/graphs/contributors)


## Computational Infrastructure for Geodynamics (CIG)

SPECFEM3D Cartesian is part of the software that is hosted by the Computational Infrastructure for Geodynamics (CIG). It is available on the CIG website [here (SPECFEM3D Cartesian)](https://geodynamics.org/resources/specfem3dcartesian).


## 2026.04.04 

* **Boundary Conditions:** Added support for **fixed** and **roller** boundary conditions within the `decompose mesh` utility and CUBIT-based mesh files.
* **Static Solver:** Implemented a static solution module using the **Conjugate Gradient (CG)** method, including support for rotational effects and other external forces.
* **Dynamic Solver:** Introduced a dynamic solution framework incorporating rotational physics (centrifugal and Coriolis forces). This includes a revised formulation for the **inverse mass matrix** to properly account for the **skew-symmetric** nature of the Coriolis matrix.
* **Examples:** Provided demonstration cases in `centri_cubit.tar.gz`. To execute the examples, follow the workflow outlined in **Steps 1–3**.