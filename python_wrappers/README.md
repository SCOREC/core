## What is _pyCore_
_pyCore_ is a python module intended to expose a minimal set of APIs from SCOREC::core so mesh adapt can be called from python. The set will include the following basic functionality:

0. Basic pcu initializations and timing functionality
1. Loading meshes and models
2. Size field computaiton
3. Call to SCOREC mesh adapt
4. Basic vtk visualization APIs

## How to build _pyCore_

### Dependencies
To build _pyCore_ the following are needed

1. [cmake](https://cmake.org/)
2. [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
3. [swig](http://www.swig.org/)
4. [SCOREC::core](git@github.com:SCOREC/core.git)

__NOTE__ If making core as a static libraries the `-fPIC` flag must be used. You can also make core as shared libraries by using the flag `-DBUILD_SHARED_LIBS=ON` option during cmake configuration. Also use the flag `-DCMAKE_SHARED_LINKER_FLAGS="-Wl,-no-as-needed` when making core with optimization on.

### Build instructions

First update `./python_wrappers/example_config.sh` to point to `SCOREC::core`'s install and then follow the steps:

1. `mkdir build`
2. `cd build`
3. `source ../example_config.sh`
4. `make`

If everything goes correctly, you will have `_pyCore.so` in your build directory. This is all you need to be able to `import` this module into python. See the example below. Note you either need to copy this file to where `PYTHONPATH` points to, or alternatively you can add the location of this file to your `PYTHONPATH`.

## How to use the _pyCore_
An example of using this module is provided in `test_pyCore.py`. You can run this code as follows

`test_pytCore.py -g ./input/cube.dmg -m ./input/cube.smb`

where `path/to/model` and `path/tp/mesh` are pointing to the location of the model and mesh files.

There is a similar example named `test_pytCore_with_simx.py` when using simmetrix:
`test_pytCore_with_simx.py -g ./input/sphere.x_t -m ./input/sphere.smb`


If Simmetrix libraries are not available update the line `-DENABLE_SIMX=ON` in `example_config.sh` to `-DENABLE_SIMX=OFF`. Note that in that case some of the functionalities of the mesh adapt will not be available or will not work as intended.

## Example models and meshes
1. `./input/cube.dmg` dmg model for the cube example
2. `./input/cube0.smb` smb (scorec format) mesh for the cube example
3. `./input/sphere.x_t` native Parasolid model the sphere example
4. `./input/sphere0.smb` smb (scorec format) mesh for the sphere example


