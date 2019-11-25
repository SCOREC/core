cmake .. \
  -DSCOREC_PREFIX="/lore/hakimm2/opt/coresimso" \ # point to the location of SCOREC/core shared libs
  -DCMAKE_INSTALL_PREFIX="../install" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DENABLE_SIMX=ON \
  -DSIM_MPI=mpich3.3 \ # has to be changed to the relevant openmpi build if using openmpi
  -DSIM_PARASOLID=ON \
  -DENABLE_PYTHON=ON
