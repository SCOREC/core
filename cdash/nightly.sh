#!/bin/bash -x
source /etc/profile
source /users/smithc11/.bash_profile

#setup lmod
export PATH=/usr/share/lmod/lmod/libexec:$PATH

#setup spack modules
unset MODULEPATH

module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
module load gcc/12.3.0-iil3lno
module load mpich/4.1.1-xpoyz4t
module load simmetrix-simmodsuite/2025.0-250108dev-llxq6sk
module load zoltan/3.83-hap4ggo
module load cmake/3.26.3-2duxfcd
module load cgns/develop-cc4dfwp

#cdash output root
d=/lore/smithc11/nightlyBuilds/
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
