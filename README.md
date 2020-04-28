# posgenpy
Python-based driver for cluster analysis using [posgen](http://apttools.sourceforge.net). These tools are to help make analysis by [posgen](http://apttools.sourceforge.net) easier.
## Installation 
...
## Overview
General description of different functions.
## Terms
**posgen** - an XML driven C++ executable which can do lots of pos-operations including generating simulated data and cluster analysis
**XML** - markup language for text file, posgen has its own specification for the structure of the file.
**cluster stats file** - the file produced by posgen which has a summary of the cluster selection results. Has cluster centres, counts of each ion (molecules not decomposed) and radius of gyration for each cluster or averaged (choice in input XML file).
## Existing MATLAB scripts
Some old MATLAB code (circa 2015) which was used by AJL previously have been uploaded, some are low level functions and others use these function to do high-level processing. There are some functions borrowed from [AtomProbeLab](http://AtomProbeLab.sourceforge.net) to let the existing functions run. These are stored in the `Matlab/` folder.
### Matlab/calculateNewClusterData.m
Low-level: Uses the [posgen](http://apttools.sourceforge.net) pos file output (`cluster.pos`, `indexed.cluster.pos`, `matrix.pos`) files to create a `.csv` file with the same format as the which IVAS generates.
Requries rangeReader, readpos (and loadMasses), ele2ionicRanges, radiusGyration, massQuant functions.
### Matlab/clusterSizeDistMultiDmax.m
High-level: Multiple random relabelling of datasets as well as a primative Dmax-sweep function, which has now been implemented in posgen. Gives an example of generating the posgen-XML file, running posgen, reading results and plotting them.
### Matlab/clusterStatsReader.m
Low-level: This reads the clusterStats.txt file produced by posgen returns the array of elements (column headings, cell array) and the decomposed counts (matrix). Uses ionStr2ions+ion2ionTable to decompose ions with some assumptions.
### Matlab/nminTable.m
Low-level: Simple function to get the number of core solute ions from the resulting cluster stats file.
### Matlab/posgenDmaxNminSweep.m
High-level: More modern script to drive posgen using the low-level `posgenOptionGen2` function and uses the sweeping capability of posgen. The function shows how to deal with the posgen output and folders, also how to organise the results.
### Matlab/posgenErosionCompSweep.m
High-level: Sweep the erosion distance and return the average cluster composition. Uses the low-level `posgenOptionGen2` function.
### Matlab/posgenOptionGen2.m
Low-level: Creates an options.xml (xmlName) file for use with posgen, includes the ability to sweep parameters.
### Matlab/posgenSweepScript2b.m
High-level: Uses Matlab/posgenErosionCompSweep to get the mean cluster composition as a function of include and erosion distance.
![Sweeping includes and erosion to get cluster composition](docs/ML_posgenSweepScript2b.PNG)
Some example results are given in `Matlab/posgenSweepScript2b_results.mat`.
### Matlab/posgenWrapperAndStats.m
High-level: Shows how to run posgen using the `posgenOptionGen2` function and make a new "IVAS like" cluster statistics file.
### Matlab/radiusGyration.m
Low-level: Helper function to calculate radius of gyration.
### Matlab/radiusGyrationMassless.m
Low-level: Helper function to calculate radius of gyration.
### Matlab/rangeReader.m
Low-level: Borrowed from [AtomProbeLab](http://AtomProbeLab.sourceforge.net) to read range files. Ideally we would use [libatomprobe](https://bitbucket.org/mycae_gmx/libatomprobe/src/default/) to read range files, but there may be another range reader already written in python we can use. Range files have no defined format and [libatomprobe](https://bitbucket.org/mycae_gmx/libatomprobe/src/default/) has good robust tools, however there are additional challenges interfacing with the C++ library.
### Matlab/solutePlotting.m
High-level: An early script which includes composition cropping as well as Nmin cropping of data.

# Authors