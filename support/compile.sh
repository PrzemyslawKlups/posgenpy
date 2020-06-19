#!/bin/sh

if [ -d libatomprobe ] ; then 
	echo "libatomprobe dir already exists. Refusing to continue";
	exit 1;
fi

if [ -d posgen ] ; then 
	echo "posgen dir already exists. Refusing to continue";
	exit 1;
fi

sudo apt-get update || { echo " Failed to update software list - internet OK?" ; exit 1 ; }

#Install build dependencies
sudo apt-get install build-essential cmake mercurial python3-dev swig || { echo "Failed to download libraries. Aborting."; exit 1 ; }

#Install actual depdendencies
sudo apt-get install libxml2-dev libgsl-dev libqhull-dev libmuparser-dev

#install depdendencies in notebooks
sudo apt-get install python3-pandas jupyter-notebook python3-matplotlib

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

make parallel || { echo "Failed to build posgen. Aborting" ; exit 1 ; } 
sudo make install || { echo "Failed to install posgen. Aborting" ; exit 1 ; }

cd ..
#=========


echo " Completed base installation. Downloading UI components"


#Download VCXSRV
#	MS announced a wayland based driver, which should eliminate the need for the below
VCXSRV=https://netcologne.dl.sourceforge.net/project/vcxsrv/vcxsrv/1.20.8.1/vcxsrv-64.1.20.8.1.installer.exe
wget "$VCXSRV" || { echo "Downloading VCXSRV failed... Skipping. IF you want a UI, you will need to do this manually" ;  VCX_DL=0 ; } ;

#This bit is to allow UI connections to work.
# You need to install an X11 server, such as vcxsrv
#Add Entry for X11 connection
echo "export DISPLAY=$(awk '/nameserver / {print $2; exit}' /etc/resolv.conf 2>/dev/null):0" >> ~/.bashrc
echo "export LIBGL_ALWAYS_INDIRECT=1" >> ~/.bashrc

source ~/.bashrc

if [ x"VCXDL" != x"0" ] ; then
	echo "Please install vcxsrv, which has been downloaded for you (5 sec)"
	sleep 5
	explorer.exe vcxsrv-64.1.20.8.1.installer.exe
else
	echo "Auto DL failed : Please manually download and install vcxsrv: https://sourceforge.net/projects/vcxsrv/ "
fi


echo "If using WSL2, and you wish to use graphical programs, Please run :"
echo "\tSet-NetFirewallProfile -DisabledInterfaceAliases 'vEthernet (WSL)' "
echo " from a windows powershell prompt. You may need to redo this for each reboot"

