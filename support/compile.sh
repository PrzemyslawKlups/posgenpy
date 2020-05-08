#!/bin/sh

if [ -d libatomprobe ] ; then 
	echo "libatomprobe dir already exists. Refusing to continue";
	exit 1;
fi

#Install build dependencies
sudo apt-get install build-essential cmake mercurial python3-dev swig || { echo "Failed to download libraries. Aborting."; exit 1 ; }

#Install actual depdendencies
sudo apt-get install libxml2-dev libgsl-dev libqhull-dev

CLONEURL=https://hg.sr.ht/~mycae/libatomprobe
hg clone $CLONEURL || { echo "Failed to obtain libraries. Aborting."; exit 1 ; }
cd libatomprobe/

#Patch to enable swig
patch -p1 < packaging/debian/patches/enable-swig.patch 

cmake . || { echo "Makefile generation failed. Aborting"; exit 1 ; }
make -j2  || { echo "library compilation failed. Aborting"; exit 1 ; }

#install library
sudo make install
