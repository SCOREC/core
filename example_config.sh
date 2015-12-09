cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
  -DENABLE_THREADS=OFF \
  -DENABLE_ZOLTAN=ON
