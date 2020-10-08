## Installation of libatomprobe for Windows Subsystem for Linux (WSL)

WSL is a light-weight virtualisation system for running linux applications
under Windows, in a "native" way. libatomprobe can be used in this way,
including with python support

To install WSL, follow these instructions, entering commands in WSL
* Install  WSL and Ubuntu App :
    Go to this URL : https://docs.microsoft.com/en-us/windows/wsl/install-win10

Run the following commands inside Ubuntu App:
* Install git:
   ```sudo apt install git```
* Inside WSL, check out this repository:
   ```git clone https://github.com/PrzemyslawKlups/posgenpy.git```
* Run the support/compile.sh script, to install all dependencies:
  first type ```$ cd posgenpy``` and then ```$ sh ./support/compile.sh```
* Run the example with :
  jupyter-notebook "Plotting real vs random clusters and more.ipynb"
 
Useful hints for using WSL: 
* paste text into the terminal using right-click (ctrl+v does not work)
* to access files from Windows drives change ```C:\...``` into ```/mnt/c/```. Remember to change all backslashes to 
    forward slashes ("\\" => "/")
* to access data from an external drive, e.g. D:
    * create a directory if needed 
        > ``` $ sudo mkdir /mnt/d ``` 
    * mount the drive
        > ```$ sudo mount -t drvfs D: /mnt/d```
    * now the data is accessable from ```/mnt/d```
    * safely unmount the drive at the end of the session 
        > ```$ sudo umount /mnt/d```

# PosgenPy
Python-based driver for cluster analysis using [posgen](http://apttools.sourceforge.net). These tools are to help make analysis by [posgen](http://apttools.sourceforge.net) easier.

## How to use PosgenPy
To make it easier and more transparent for you to analyse your clusters, we created some Jupyter Notebook files in which you can see our examples, tweak them according to your preference or insert path files to your `pos` and `rrng` files and use it directly on your data.
We recommend starting from `Quick_Parameter_Check.ipynb`. In this notebook, you will be able to run an analysis on your data by substituting file paths and adjusting parameters if required.

## Overview
General description of different functions.

**write_xml_with_relabelling** Function to write an XML file that is read by `posgen` to perform max-sep clustering on a `pos` file. It includes relabelling function which helps with estimating the number of real clusters in the data.


## Terms
**posgen** - an XML driven C++ executable which can do lots of pos-operations including generating simulated data and cluster analysis

**XML** - markup language for text file, posgen has its own specification for the structure of the file.

**cluster stats file** - the file produced by posgen which has a summary of the cluster selection results. Has cluster centres, counts of each ion (molecules not decomposed) and radius of gyration for each cluster or averaged (choice in input XML file).

# Authors
PosgenPy was created by Dan Hayley, Andy London, James Famelton, Hazel Gardner, and Przemek Klups.
