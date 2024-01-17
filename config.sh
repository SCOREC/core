cmake -S . -B /lore/mccalf/build \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DIS_TESTING=on \
    -DBUILD_EXES=on \
    -DMESHES=/users/mccalf/pcu-update/pumi-meshes
