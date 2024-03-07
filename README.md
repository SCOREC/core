# SCOREC Core #


The SCOREC Core is a set of C/C++ libraries for unstructured mesh
simulations on supercomputers.

For more information, start at our
[wiki page](https://github.com/SCOREC/core/wiki)

### What is in this repository? ###

* PUMI: parallel unstructured mesh infrastructure API
  [User's Guide](http://scorec.rpi.edu/~seol/PUMI.pdf)
* PCU: Communication and file IO built on MPI 
* APF: Abstract definition of meshes, fields, and related operations
* GMI: Common interface for geometric modeling kernels
* MDS: Compact but flexible array-based mesh data structure
* PARMA: Scalable partitioning and load balancing procedures
* SPR: Superconvergent Patch Recovery error estimator
* MA: Anisotropic mixed mesh adaptation and solution transfer
* SAM: Sizing anisotropic meshes
* STK: Conversion from APF meshes to Sandia's STK meshes
* ZOLTAN: Interface to run Sandia's Zoltan code on APF meshes
* PHASTA: Tools and file formats related to the PHASTA fluid solver
* MTH: Math containers and routines
* CRV: Support for curved meshes with Bezier Shapes
* PYCORE: Python Wrappers (see python_wrappers/README.md for build instructions)
* REE: Residual based implicit error estimator

### How do I get set up? ###

* Dependencies: CMake for compiling and MPI for running
* Configuration: Typical CMake configure and build.
  The `example_config.sh` shows common options to select,
  use a front-end like `ccmake` to see a full list of options
* Tests: the test/ subdirectory has tests and standalone
  tools that can be compiled by explicitly listing them as targets
  to `make`.
* Users: `make install` places libraries and headers in
  a specified prefix, application code can use these
  in their own compilation process.
  We also install pkg-config files for all libraries.

### Contribution guidelines ###

* Don't break the build
* See the `STYLE` file
* If in doubt, make a branch
* Run the ctest suite
* Don't try to force push to `master` or `develop`; it is disabled

### Who do I talk to? ###

* If you have a usage question or have found a bug please post an [issue](https://github.com/SCOREC/core/issues).
* Otherwise, email [Cameron Smith](https://www.scorec.rpi.edu/~cwsmith) and pumi@scorec.rpi.edu.

### Citing PUMI

If you use these tools, please cite the following paper:

Daniel A. Ibanez, E. Seegyoung Seol, Cameron W. Smith, and Mark S. Shephard. 2016. PUMI: Parallel Unstructured Mesh Infrastructure. ACM Trans. Math. Softw. 42, 3, Article 17 (May 2016), 28 pages. DOI: https://doi.org/10.1145/2814935

We would be happy to provide feedback on journal submissions using PUMI prior to publication.
