#!/bin/bash

## updates list of libraries & install Boost
sudo apt-get update
sudo apt-get install libboost-all-dev

## install freesasa and its prerequisities autoconf and git
sudo apt-get install git
sudo apt-get install autoconf
git clone https://github.com/mittinatten/freesasa.git
cd freesasa
autoreconf -i
./configure --disable-xml --disable-json
make
sudo make install
cd ..

## compile inspire
cd src
make
## If you want to install all tools, use the following instead
# make all
sudo make install
