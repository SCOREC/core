cmake .. \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_C_FLAGS="-Wall -g -O2" \
  -DCMAKE_CXX_FLAGS="-Wall -g -O2" \
  -DENABLE_THREADS=ON \
  -DSIM_MPI="openmpi1.6.5" \
