# # # #
# The following line is for capstone compilation
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/cwsmith/projects/pettt/capstone/rhel7/v900rc2/lib
# # # #
export MESHES=/lore/woodra/pumi-meshes

flags="-g -O0"
cmake .. \
  -DCMAKE_BUILD_TYPE="Debug" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_EXE_LINKER_FLAGS="-lpthread ${flags}" \
  -DSCOREC_CXX_OPTIMIZE=ON \
  -DBUILD_EXES=ON \
  -DMDS_ID_TYPE="int" \
  -DENABLE_CAPSTONE=ON \
  -DIS_TESTING=ON \
  -DMESHES="$MESHES" \
  -DENABLE_ZOLTAN=ON \
  -DPCU_COMPRESS=ON \
  -DENABLE_FPP=ON


