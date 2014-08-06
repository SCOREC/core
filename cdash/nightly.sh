#!/bin/bash

source /usr/local/etc/bash_profile
module load cmake/latest
module load openmpi/1.6.5-ib
module load parmetis/OpenMPI-1.6.5/4.0.2
module load zoltan/OpenMPI-1.6.5/3.8
module load simmetrix/simModSuite
module load netcdf

export PATH=/lore/dibanez/bin:$PATH

cd /lore/dibanez/cdash

svn up meshes/
ctest -VV -D Nightly -S /lore/dibanez/core/cdash/nightly.cmake &> log

rm -rf build/
