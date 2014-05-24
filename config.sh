cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_FLAGS="-Wall -g -O2" \
  -DCMAKE_CXX_FLAGS="-Wall -g -O2" \
  -DENABLE_THREADS=ON \
