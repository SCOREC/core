PETSC_DIR=/fasttmp/seol/petsc-3.5.4
PETSC_ARCH=real-mpich3-gcc4.4.5
ZOLTAN_DIR=/fasttmp/seol/mpich3-gcc4.4.5-install
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/mpich3/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/mpich3/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O2" \
  -DCMAKE_CXX_FLAGS=" -g -O2" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX=$ZOLTAN_DIR \
  -DENABLE_OMEGA_H=ON \
  -DCMAKE_BUILD_TYPE=Debug

