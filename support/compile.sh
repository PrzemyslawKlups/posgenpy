#!/bin/sh

if [ -d libatomprobe ] ; then 
	echo "libatomprobe dir already exists. Refusing to continue";
	exit 1;
fi

if [ -d posgen ] ; then 
	echo "posgen dir already exists. Refusing to continue";
	exit 1;
fi

#Install build dependencies
sudo apt-get install build-essential cmake mercurial python3-dev swig || { echo "Failed to download libraries. Aborting."; exit 1 ; }

#Install actual depdendencies
sudo apt-get install libxml2-dev libgsl-dev libqhull-dev libmuparser-dev

#install depdendencies in notebooks
sudo apt-get install python3-pandas jupyter-notebook

#Download and compile libatomprobe
#=====
CLONEURL=https://hg.sr.ht/~mycae/libatomprobe
hg clone $CLONEURL || { echo "Failed to obtain libraries. Aborting."; exit 1 ; }
cd libatomprobe/

#Patch to enable swig
patch -p1 < packaging/debian/patches/enable-swig.patch 

cmake . || { echo "Makefile generation failed. Aborting"; exit 1 ; }
make  || { echo "library compilation failed. Aborting"; exit 1 ; }

#install library
sudo make install
#Run linker configurer
sudo ldconfig

cd ..
#=========

#Download and compile posgen
#==========
POSGEN_CLONE=http://hg.code.sf.net/p/apttools/posgen/code 
hg clone $POSGEN_CLONE posgen || { echo "Failed to download posgen source code. Aborting."; exit 1 ; }

cd posgen

make || { echo "Failed to build posgen. Aborting" ; exit 1 ; } 
sudo make install || {echo "Failed to install posgen. Aborting" ; exit 1 ; }

cd ..
#=========


echo " Completed installation."
