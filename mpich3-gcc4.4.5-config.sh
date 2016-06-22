PETSC_DIR=/fasttmp/seol/petsc-3.5.4
PETSC_ARCH=real-mpich3-gcc4.4.5

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/mpich3/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/mpich3/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O2" \
  -DCMAKE_CXX_FLAGS=" -g -O2" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="/fasttmp/seol/mpich3-gcc4.4.5-install/include" \
  -DZOLTAN_LIBRARY="/fasttmp/seol/mpich3-gcc4.4.5-install/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PETSC_DIR/$PETSC_ARCH/include" \
  -DPARMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PETSC_DIR/$PETSC_ARCH/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX="/fasttmp/seol/mpich3-gcc4.4.5-install" \
  -DENABLE_OMEGA_H=ON \
  -DCMAKE_BUILD_TYPE=Debug

