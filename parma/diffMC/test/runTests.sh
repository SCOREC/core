#!/bin/bash -e
mpirun=/usr/local/openmpi/latest/bin/mpirun 
declare -a prio=("-v 1" "-r 1" "-r 2 -v 1" "-r 1 -v 2")
for p in "${prio[@]}"; do 
   pNoSpace=`echo $p | tr -d ' '`
   logfile="logAfosr${pNoSpace}"
   echo $logfile
   set -x
   $mpirun -np 16 ./diffMcDriver -m /fasttmp/cwsmith/meshes/svnMeshes/airFoilAfosr/16/afosr_60e3tets_16p.sms -t 0 -i 20 $p &> $logfile  
   set +x
   grep 'STATUS ParMA elapsed time' $logfile
   grep 'number of disconnected components' $logfile
   grep '^imbalance' $logfile
done
#$mpirun -np 2 ./diffMcDriver -m geom_uni_2_.sms -t 3 -i 10 &> logt3 &
