# Disco4est (d4est)

A scalable hp-adaptive discontinuous Galerkin solver for coupled non-linear elliptic partial differential equations on curved multi-block meshes.

### Prerequisites

The dependencies are OpenBLAS, HDF5, PETSc and p4est. You will need at least OpenBLAS and HDF5. If you do not have PETSc or p4est, these come bundled, so they can be installed automatically.

### Installing

1) git clone --recursive https://github.com/trevor-vincent/d4est d4est

2) mkdir build

3) Open d4est/Support/machines.cmake and edit the LAPTOP if in machines.cmake to fit your needs
If you have already installed p4est, zlib and petsc, please look at LAPTOP_WITH_PREFIX instead

4) cd build && cmake ../d4est -DLAPTOP -DCMAKE_BUILD_TYPE=DEBUG && make -j4

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
