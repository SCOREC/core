#!/bin/bash -x

# load environment variables
source /usr/local/etc/bash_profile
module load cmake/latest
module load mpich3/3.1.2-thread-multiple
module load parmetis/mpich3.1.2/4.0.3
module load zoltan/mpich3.1.2/3.8
module load simmetrix/simModSuite
module load netcdf
module load llvm/latest
module load git

#cdash output root
cd /lore/dibanez/cdash

#run nightly.cmake script
ctest -VV -D Nightly -S /lore/dibanez/core/cdash/nightly.cmake &> cmake_log.txt

#core repository checked out by nightly.cmake
cd /lore/dibanez/cdash/build/core
#build the Doxygen html documentation
make doc
#remove the old web documentation
rm -rf /net/web/public/dibanez/core
#replace it with the generated one
cp -r doc/html /net/web/public/dibanez/core

#core repository checked out by nightly.cmake
cd /lore/dibanez/cdash/build/core
#clean the build of object files
make clean
#run Coverity static analysis on the build
export PATH=$PATH:/lore/dibanez/cov-analysis-linux64-7.7.0/bin
cov-build --dir cov-int make -j 4
#pack up the tarball of results
tar czvf pumi.tgz cov-int
#send it to Coverity Scan site using curl
curl --form token=faMZVmxTjByhNoJyb_4wDw \
  --form email=dan.a.ibanez@gmail.com \
  --form file=@$PWD/pumi.tgz \
  --form version="1.0.1" \
  --form description="Automated" \
  https://scan.coverity.com/builds?project=SCOREC%2Fcore

#remove compilation directories created by nightly.cmake
cd /lore/dibanez/cdash
rm -rf build/

cd /lore/dibanez
find meshes/phasta -name "*procs_case" | xargs rm -rf
find meshes/phasta -name "out_mesh" | xargs rm -rf
