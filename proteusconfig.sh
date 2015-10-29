cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-O2 -g -Wall -fPIC" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall -fPIC" \
  -DENABLE_ZOLTAN=ON \
  -DSIM_MPI="mpich3" \
  -DENABLE_THREADS=ON
