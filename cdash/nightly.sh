#!/bin/bash -x
source /etc/profile
source /users/cwsmith/.bash_profile

#setup lmod
export PATH=/usr/share/lmod/lmod/libexec:$PATH

#setup spack modules
unset MODULEPATH

module unuse /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core
module use /opt/scorec/modules
module use /opt/scorec/spack/v0154_2/lmod/linux-rhel7-x86_64/Core
module load gcc/10.1.0
module load mpich/3.3.2
module load simmetrix-simmodsuite/17.0-220516
module load zoltan/3.83-int32
module load cmake/3.20.0

#cdash output root
d=/lore/cwsmith/nightlyBuilds/
cd $d
#remove compilation directories created by previous nightly.cmake runs
[ -d build ] && rm -rf build/

touch $d/startedCoreNightly
#run nightly.cmake script
ctest -V --script $d/nightly.cmake
touch $d/doneCoreNightly

#create doxygen docs
cd build/master
make doc
cp -r doc/html /net/web/scorec/scorec-web/htdocs/pumi/doxygen
