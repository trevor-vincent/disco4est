[![license](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://github.com/trevor-vincent/d4est/blob/master/LICENSE.md)
[![Testing (Linux)](https://github.com/trevor-vincent/disco4est/actions/workflows/tests_linux.yml/badge.svg)](https://github.com/trevor-vincent/disco4est/actions/workflows/tests_linux.yml)
[![arXiv](https://img.shields.io/badge/arXiv-1907.01572-green.svg?style=flat)](https://arxiv.org/abs/1907.01572)
<p align="center">
<img src="https://github.com/trevor-vincent/d4est/blob/master/support/Logo/logo.png" width="250">
</p>



disco4est (abbrev. d4est) is a scalable hp-adaptive discontinuous Galerkin solver for non-linear elliptic partial differential equations on curved multi-block meshes.

## Solver Details

<img src="https://github.com/trevor-vincent/disco4est/blob/master/2pun.png" width="400px" align="right">

- The **disco**ntinuous Galerkin method is used to discretize the non-linear or linear elliptic PDEs.

- We use an in-house multigrid algorithm to precondition PETSc solvers. 

- The multigrid algorithm can handle different smoothers (Chebyshev, Additive Schwarz... ) and different bottom solvers (p-multigrid, Krylov, ...). 

- The multigrid algorithm uses surrogate intermediary meshes to keep the cells 2:1 balanced.

- Different topologies and parallel mesh refinement is handled using the multi-block tree code p4est. 

- Scaling of the multigrid solver has been tested up to O(10,000) cores.

## Prerequisites

The dependencies are MPI, OpenBLAS, HDF5, PETSc, zlib, zlog and p4est. You will need at least MPI. If you do not have zlib, zlog, PETSc p4est, HDF5 or OpenBLAS, then these come bundled, so they will be installed automatically to the local build folder unless you specify otherwise.

> On Mac OS X, you may simply install most dependencies through Homebrew:
>
> ```bash
> brew install openblas hdf5 gperftools zlib
> ```
>
> Just remember that the OS X system `gcc` is just a `clang` wrapper, so make sure to set the `CC` and `CXX` environment variables to e.g. `gcc-7` and `g++-7`, with `gcc` (along with its bundled `gfortran`) installed through Homebrew. Also make sure your MPI installation wraps the same `gcc`, e.g. through `brew install mpich --cc=gcc-7`.
>
> Note that due to a bug in the `syslog` headers in OS X 10.13 (High Sierra) you should provide the flag `mmacosx-version-min`, such as in `export CC="gcc-7 -mmacosx-version-min=10.12"` to compile `zlog`.

## Installation

1) `git clone --recursive https://github.com/trevor-vincent/d4est && cd d4est`

2) Duplicate one of the `Support/CMake/machine.cmake.*.example` build configuration files to `Support/CMake/machine.cmake` and modify it to fit the local machine.
  - Adjust the `CMAKE_BUILD_TYPE` variable to `Release` for (much faster) production builds.
  - If you choose to use the bundled libraries, then these will be located in the build directory and will be compiled only once, unless you change cmake commandline options or delete the build directory.
  - The cmake script searches for the MPI build on your system, so this never needs to be specified.
  - There is a preconfigured example file for the _Minerva_ cluster. On Minerva, also `source` `Support/Scripts/build.env.minerva.example` to configure the build environment, or duplicate and modify the file first.

3) Compile:

  ```bash
  mkdir build && cd build
  cmake ..
  make -j4 # or adjust number of cores
  ```

4) Run an example:

  ```
  cd ConstantDensityStar
  mpirun -np 4 path/to/build/ConstantDensityStar/constant_density_star_driver
  ```

  For the example to complete in O(10s) you may want to consult the `options.input` file to set `amr.num_of_amr_steps=1`.
  
  On a cluster, refer to the `Support/Scripts/submit.*.*.example` files for guidance.


## Options documentation

Configuration options are loaded from an `options.input` file residing in the same directory as the problem executable by default, or from the file supplied as the first positional command line argument.

- `[problem]`: Problem-specific options.
- `[initial_mesh]`
  - `min_quadrants`
  - `min_level`: Number of initial mesh refinements.
  - `fill_uniform`
  - `regionX_deg`: Initial polynomial degree in region X.
  - `regionX_deg_quad_inc`: Initial quadrature polynomial degree increment in region X.
- `[amr]`
  - `scheme` can be one of the following AMR refinement schemes:
    - `uniform_h`: Divides each element into 4 sub-elements in 2 dimensions, or into 8 sub-elements in 3 dimension.
    - `uniform_p`: Increases the polynomial order on each element by 1.
    - `smooth_pred`: hp-amr scheme that uses estimator performance to determine h or p refinement
  - `num_of_amr_steps`
  - `max_degree`

## Licenses

This project is licensed under the GPLv2 License - see the [LICENSE.md](LICENSE.md) file for details
</a><br />The logo is licensed under a <a rel="license" href="https://creativecommons.org/licenses/by/3.0/">Creative Commons  License</a>.

The disco-ball in the logo is a modified version of a design by Daouna Jeong from the Noun Project and the tree in the logo is a modified version of a design created by Musket from the Noun Project.
