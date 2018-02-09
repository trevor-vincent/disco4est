[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/trevor-vincent/d4est/blob/master/LICENSE.md)

# disco4est (d4est)

A scalable hp-adaptive discontinuous Galerkin solver for coupled non-linear elliptic partial differential equations on curved multi-block meshes.

## Prerequisites

The dependencies are OpenMPI, OpenBLAS, HDF5, PETSc, zlib and p4est. You will need at least OpenMPI, OpenBLAS and HDF5. If you do not have zlib, PETSc or p4est, these come bundled, so they will be installed automatically to the local build folder unless you specify otherwise.

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

2) Duplicate `Support/machine.cmake.example` to `Support/machine.cmake` and edit the paths to the dependencies.
  - If you have already installed p4est, zlib and petsc, also give their respective paths.
  - Adjust the `CMAKE_BUILD_TYPE` variable to `Release` for (much faster) production builds.
  - If you choose to use the bundled PETSc, p4est and zlib, then these will be located in the build directory and will be compiled only once, unless you change cmake commandline options or delete the build directory.
  - The cmake script searches for the MPI build on your system, so this never needs to be specified. The cmake script does not search for any of the other dependencies (this is a TODO).

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

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
