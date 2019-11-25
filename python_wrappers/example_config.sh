# DSCOREC_PREFIX should point to the location of SCOREC/core shared libs
# DSIM_MPI should be set to relevant mpi build
cmake .. \
  -DSCOREC_PREFIX="/lore/hakimm2/opt/coresimso" \
  -DCMAKE_INSTALL_PREFIX="../install" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DENABLE_SIMX=ON \
  -DSIM_MPI=mpich3.3 \
  -DSIM_PARASOLID=ON \
  -DENABLE_PYTHON=ON
