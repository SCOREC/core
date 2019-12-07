## What is _pyCore_
_pyCore_ is a python module intended to expose a minimal set of APIs from SCOREC::core so mesh adapt can be called from python. The set will include the following basic functionality:

0. Basic pcu initializations and timing functionality
1. Loading meshes and models
2. Size field computation
3. Call to SCOREC mesh adapt
4. Basic vtk visualization APIs

## How to build _pyCore_

### Dependencies
To build _pyCore_ the following are needed

1. [cmake](https://cmake.org/)
2. [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
3. [swig](http://www.swig.org/)
4. [SCOREC::core](git@github.com:SCOREC/core.git)

__NOTE__ If making core as a static libraries the `-fPIC` flag must be used. You can also make core as shared libraries by using the flag `-DBUILD_SHARED_LIBS=ON` option during cmake configuration.

### Build instructions

First update `core/example_config_with_python_interface.sh`. The main options that you need to set are `-DPUMI_PYTHON_INTERFACE` to `ON`. To build with simmetrix tool also `-DENABLE_SIMMETRIX` to `ON`, `-DSIM_MPI=mpi_version_used_by_sim_libs`, and `-DSIM_PARASOLID` or `-DSIM_ACICS` to `ON` depending of the CAD kernel you are using. Then follow the instructions for the building and installing SCOREC::core.

If everything goes correctly, you will have `_pyCore.so` and `pyCore.py` (and also `simhelper.so` if you build with simmetrix enables) in your install directory specified by `-DCMAKE_INSTALL_PREFIX`. This is all you need to be able to `import` this module into python. See the example below. Note that will need to add the install location of these files to your `PYTHONPATH`.

## How to use the _pyCore_
An example of using this module is provided in `test_pyCore.py`. You can run this code as follows

`test_pytCore.py -g ./input/cube.dmg -m ./input/cube.smb`

where `path/to/model` and `path/tp/mesh` are pointing to the location of the model and mesh files.

There is a similar example named `test_pytCore_with_simx.py` when using simmetrix:
`test_pytCore_with_simx.py -g ./input/sphere.x_t -m ./input/sphere.smb`


## Example models and meshes
1. `./input/cube.dmg` dmg model for the cube example
2. `./input/cube0.smb` smb (scorec format) mesh for the cube example
3. `./input/sphere.x_t` native Parasolid model the sphere example
4. `./input/sphere0.smb` smb (scorec format) mesh for the sphere example

