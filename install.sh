#!/usr/bin/bash

#make rabbitFX library
cd RabbitFX &&
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install && 
cd ../../ &&

#compile RabbitKSSD
mkdir -p build
cd build
cmake ..
make -j8 && make install
cd ..
