#!/bin/bash -ex
module load simmetrix/simModSuite/9.0-140403

g++ \
-g \
$PWD/cadToSim.cc \
-o cadToSim \
-I $SIM_INCLUDE_DIR \
-L $SIM_LIB_DIR \
-lSimMeshTools \
-lSimParasolid260 \
-lSimAcis240 \
-lSimModel \
$SIM_LIB_DIR/psKrnl/libpskernel.so \
$SIM_LIB_DIR/acisKrnl/libSpaACIS.so 
