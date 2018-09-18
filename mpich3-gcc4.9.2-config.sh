PARMETIS_DIR=/lore/esyoon/openlib/petsc-3.7.5-mpich3-gcc-4.9.2
ZOLTAN_DIR=/lore/seol/mpich3-gcc4.9.2-install
PREFIX=/lore/seol/mpich3-gcc4.9.2-install

cmake .. \
  -DCMAKE_C_COMPILER="/usr/local/mpich3/latest/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="/usr/local/mpich3/latest/bin/mpicxx" \
  -DCMAKE_C_FLAGS=" -g -O0 -DDEBUG -Wall -Wextra -Werror" \
  -DCMAKE_CXX_FLAGS=" -g -O0 -DDEBUG -Wall -Wextra -Werror" \
  -DENABLE_ZOLTAN=ON \
  -DZOLTAN_INCLUDE_DIR="$ZOLTAN_DIR/include" \
  -DZOLTAN_LIBRARY="$ZOLTAN_DIR/lib/libzoltan.a" \
  -DPARMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DPARMETIS_LIBRARY="$PARMETIS_DIR/lib/libparmetis.so" \
  -DMETIS_INCLUDE_DIR="$PARMETIS_DIR/include" \
  -DMETIS_LIBRARY="$PARMETIS_DIR/lib/libmetis.so" \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DSCOREC_CXX_WARNINGS=ON \
  -DBUILD_EXES=ON \
  -DIS_TESTING=ON \
  -DMESHES=/fasttmp/seol/scorec/meshes \
  -DMPIRUN=/usr/local/mpich3/latest/bin/mpirun \
  -DCMAKE_BUILD_TYPE=Debug

