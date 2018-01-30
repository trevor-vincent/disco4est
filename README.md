# disco4est (d4est)

A scalable hp-adaptive discontinuous Galerkin solver for coupled non-linear elliptic partial differential equations on curved multi-block meshes.

### Prerequisites

The dependencies are OpenMPI, OpenBLAS, HDF5, PETSc, zlib and p4est. You will need at least OpenMPI, OpenBLAS and HDF5. If you do not have zlib, PETSc or p4est, these come bundled, so they will be installed automatically to the local build folder unless you specify otherwise. 



### Installing

1) git clone --recursive https://github.com/trevor-vincent/d4est

2) mkdir build

3) Open d4est/Support/cmake/machines.cmake and edit the LAPTOP section in machines.cmake (or make a new section) to fit your needs. If you have already installed p4est, zlib and petsc on your system please look at the LAPTOP_WITH_PREFIX section in machines.cmake instead of LAPTOP.

4) cd build && cmake ../d4est -DLAPTOP=1 -DCMAKE_BUILD_TYPE=DEBUG && make -j4

Note 1): -DCMAKE_BUILD_TYPE=Release for production builds (much faster) and make -j{cores_on_your_machine} if you have more than 4 cores

Note 2): If you choose to use the bundled PETSc, p4est, zlib, then these will be located in the build directory and will be compiled only once, unless you change cmake commandline options or delete the build directory. 

Note 3): The cmake script searches for the MPI build on your system, so this never needs to be specified. The cmake script does not search for any of the other dependencies (this is a TODO).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
