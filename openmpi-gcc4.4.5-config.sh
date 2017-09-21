PETSC_DIR=/lore/seol/petsc-3.5.4
PETSC_ARCH=real-openmpi1.6.5
PARMETIS_DIR=$PETSC_DIR/$PETSC_ARCH
PREFIX=/lore/seol/openmpi-gcc4.4.5-install
ZOLTAN_DIR=$PREFIX
cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/openmpi/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/openmpi/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O2 -Wall -Wextra -Werror" \
  -DCMAKE_CXX_FLAGS=" -g -O2 -Wall -Wextra -Werror" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.a" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.a" \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DSCOREC_CXX_WARNINGS=ON \
  -DBUILD_EXES=ON \
  -DIS_TESTING=ON \
  -DMESHES=/fasttmp/seol/scorec/meshes \
  -DMPIRUN=/usr/local/openmpi/latest/bin/mpirun \
  -DCMAKE_BUILD_TYPE=Debug

