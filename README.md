[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/trevor-vincent/d4est/blob/master/LICENSE.md)
[![Build Status](https://travis-ci.org/trevor-vincent/d4est.svg?branch=master)](https://travis-ci.org/trevor-vincent/d4est)
<p align="center">
<img src="http://cita.utoronto.ca/~tvincent/logo.png" width="250">
</p>
disco4est (d4est) is a scalable hp-adaptive discontinuous Galerkin solver for coupled non-linear elliptic partial differential equations on curved multi-block meshes.

## Prerequisites

The dependencies are MPI, OpenBLAS, HDF5, PETSc, zlib and p4est. You will need at least MPI. If you do not have zlib, PETSc p4est, HDF5 or OpenBLAS, then these come bundled, so they will be installed automatically to the local build folder unless you specify otherwise.

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

2) Duplicate `Support/machine.cmake.prefix.example` or `Support/machine.cmake.bundled.example` to `Support/machine.cmake` and edit the paths to the dependencies.
  - If you have already installed p4est, zlib and petsc, also give their respective paths.
  - Adjust the `CMAKE_BUILD_TYPE` variable to `Release` for (much faster) production builds.
  - If you choose to use the bundled libraries, then these will be located in the build directory and will be compiled only once, unless you change cmake commandline options or delete the build directory.
  - The cmake script searches for the MPI build on your system, so this never needs to be specified.

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


## Options documentation

Configuration options are loaded from an `options.input` file residing in the same directory as the problem executable by default, or from the file supplied as the first positional command line argument.

- `[problem]`: Problem-specific options.
- `[initial_mesh]`
  - `min_quadrants`
  - `min_level`: Number of initial mesh refinements.
  - `fill_uniform`
  - `regionX_deg`: Initial polynomial degree in region X.
  - `region0_deg_quad_inc`
- `[amr]`
  - `scheme` can be one of the following AMR refinement schemes:
    - `uniform_h`: Divides each element into 4 sub-elements in 2 dimensions, or into 8 sub-elements in 3 dimension.
    - `uniform_p`: Increases the polynomial order on each element by 1.
  - `num_of_amr_steps`
  - `max_degree`


## Licenses

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
</a><br />The logo is licensed under a <a rel="license" href="https://creativecommons.org/licenses/by/3.0/">Creative Commons  License</a>.

The disco-ball in the logo is a modified version of a design by Daouna Jeong from the Noun Project and the tree in the logo is a modified version of a design created by Musket from the Noun Project.
