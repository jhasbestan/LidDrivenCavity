#! /usr/bin/bash

module use /mnt/share/intel/oneapi/2022.2/modulefiles
module load icc/2022.1.0

export KMP_AFFINITY=granularity=fine,compact

export CXX=icpc
export CC=icc

module list 

cd build && rm CMakeCache.txt


cmake ../ -DCMAKE_BUILD_TYPE=Release -DZMM_HIGH=true  && make -j 8  && make install

cd ../bin

./incomp_avx512 40000 10

cd ../build && rm CMakeCache.txt

cmake ../ -DCMAKE_BUILD_TYPE=Release -DAVX2=true  && make -j 8 && make install 

cd ../bin

./incomp_avx2 40000 10

cd ../

find . -iname *.optrpt

