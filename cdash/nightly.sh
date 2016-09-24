#!/bin/bash -x

# load environment variables
source /usr/local/etc/bash_profile
module load cmake/latest
module load mpich3/3.1.2-thread-multiple
module load parmetis/mpich3.1.2/4.0.3
module load zoltan/mpich3.1.2/3.8
module load simmetrix/simModSuite/11.0-160922pre
module load netcdf
module load gcc/4.9.2
module load git
module load valgrind/3.8.1

#cdash output root
cd /lore/dibanez/cdash
#remove compilation directories created by previous nightly.cmake runs
rm -rf build/

#run nightly.cmake script
ctest --output-on-failure --script /lore/dibanez/src/core/cdash/nightly.cmake &> cmake_log.txt
cp cmake_log.txt /net/web/public/dibanez/nightly_cmake_log.txt

if [ -d "/lore/dibanez/cdash/build/master" ]; then
  #core repository checked out by nightly.cmake
  cd /lore/dibanez/cdash/build/master
  #build the Doxygen html documentation
  make doc
  if [ -d "$PWD/doc/html" ]; then
    #remove the old web documentation
    rm -rf /net/web/public/dibanez/core
    #replace it with the generated one
    cp -r doc/html /net/web/public/dibanez/core
  fi
fi

#core repository checked out by nightly.cmake
cd /lore/dibanez/cdash/build/master
#clean the build of object files
make clean
#run Coverity static analysis on the build
export PATH=$PATH:/lore/dibanez/cov-analysis-linux64-7.7.0.4/bin
cov-build --dir cov-int make -j 4
#pack up the tarball of results
tar czvf pumi.tgz cov-int
#cleanup the Chef test output
cd /lore/dibanez
find meshes/phasta -name "*procs_case" | xargs rm -rf
find meshes/phasta -name "out_mesh" | xargs rm -rf
