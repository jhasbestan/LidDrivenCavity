#! /usr/bin/bash

module use /mnt/share/intel/oneapi/2022.2/modulefiles
module load icc/2022.1.0


export CXX=icpc
export CC=icc

module list 

cd build

cmake ../ -DCMAKE_BUILD_TYPE=Release && make install 

cd ../bin


export KMP_AFFINITY=granularity=fine,compact

./incomp 200




