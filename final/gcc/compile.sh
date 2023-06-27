#! /usr/bin/bash

module use /mnt/share/gcc-13.1/modulefiles
module load gcc-13.1 

export CXX=g++
export CC=gcc

echo "====================================================================================" | tee -a compile.log
printf "\n \n"  | tee -a compile.log

build=buildavx512
if [[ ! -e $build ]]; then
    echo "making build avx512 directory-------------------" | tee -a compile.log
    mkdir $build
else
    echo "cleaning up the build directory-------------------" | tee -a compile.log
    rm -rf $build/*
fi

cd $build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DZMM_HIGH=true  && make -j8

cd ..
build=buildavx2 
if [[ ! -e $build ]]; then
    echo "making build avx2 directory-------------------" | tee -a compile.log
    mkdir $build
else
    echo "cleaning up the build directory-------------------" | tee -a compile.log
    rm -rf $build/*
fi

cd $build

cmake ../ -DCMAKE_BUILD_TYPE=Release -DAVX2=true  && make -j8

