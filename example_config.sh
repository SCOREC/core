cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DSIM_MPI="mpich3.1.2" \
  -DENABLE_ZOLTAN=ON
