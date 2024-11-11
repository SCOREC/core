#!/bin/bash -e

branch=$1

cd $SCRATCH/globus-compute/core-test

# core
rm build-core -rf
rm core -rf
git clone https://github.com/SCOREC/core.git
cd core && git checkout $branch && git clone https://github.com/SCOREC/pumi-meshes.git && cd -
cmake -S core -B build-core \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DENABLE_ZOLTAN=OFF \
  -DMPIRUN=/usr/bin/srun \
  -DMPIRUN_PROCFLAG="--ntasks" \
  -DIS_TESTING=True
cmake --build build-core -j 24