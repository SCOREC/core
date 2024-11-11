#!/bin/bash

name=core

cd $SCRATCH/globus-compute/$name-test

cd build-$name
salloc --time 00:20:00 --constrain=gpu --qos=interactive --nodes=1 --ntasks-per-node=40 --cpus-per-task=1 --gpus=1 --account=m4564 ctest
cat $PWD/Testing/Temporary/LastTest.log