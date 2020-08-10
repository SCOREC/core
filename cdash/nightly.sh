#!/bin/bash -x
source /etc/profile
source /users/cwsmith/.bash_profile

#setup lmod
export PATH=/usr/share/lmod/lmod/libexec:$PATH

#setup spack modules
unset MODULEPATH
module use /opt/scorec/spack/lmod/linux-rhel7-x86_64/Core/

module load gcc/7.3.0-bt47fwr
module load cmake/3.12.1-wfk2b7e
module load mpich
module load zoltan/3.83-int32-gx3prjr
module load simmetrix-simmodsuite/14.0-190928-snlypqg

#cdash output root
d=/lore/cwsmith/nightlyBuilds/
cd $d
#remove compilation directories created by previous nightly.cmake runs
[ -d build ] && rm -rf build/

touch $d/startedCoreNightly
#run nightly.cmake script
ctest -V --script $d/repos/core/cdash/nightly.cmake
touch $d/doneCoreNightly

#create doxygen docs
cd build/master
make doc
cp -r doc/html/* /net/web/scorec/scorec-web/htdocs/pumi/doxygen/.
