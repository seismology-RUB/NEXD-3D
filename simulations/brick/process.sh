#!/bin/bash
#

###################################################

echo "running example: `date`"
currentdir=`pwd`
sourcedir=~/git/dg3d/

# sets up directory structure in current directory
echo
echo "   setting up example..."
echo

mkdir -p bin
mkdir -p out

rm -rf out/*
rm -f run.csh.*

# compiles executables in root directory
cd $sourcedir
#make clean
make all > tmp.log
cd $currentdir

# copy program
cd bin/
cp ${sourcedir}bin/mesher .
cp ${sourcedir}bin/solver .
cp ${sourcedir}bin/movie .
cd ../

# mesh
./bin/mesher

# stores setup
cp data/parfile out/
cp data/stations out/
cp data/source out/



