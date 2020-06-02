%Scripts, functions and XML templates written by James Famelton for processing of APT data. 
Sections are written or adapted from the work of Andrew J. London and Daniel Haley and are creditted in the code were approriate

This software is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

REQUIREMENTS AND INSTALL

Some functions and scripts require the use of the command line progam posgen. Avaialbe from http://apttools.sourceforge.net/
The use of posgen requires the Matlab working directory to be the folder which contains posgen.exe and posscript.dtd.
The aim is for all these scripts to be directly copyable into that directory. 
Some functions will also look for XML templates to feed to posgen from the same directory. 

XML TEMPLATES

Owing to a lack of skill and not wanting to waste time on authors part, the Matlab code does not write it's own XML file for posgen.
Instead it edits files that already exist and either saves a new copy or overwrites. There are four such functions;
bulkcounts
cluster_make_pos
ClusterSweep_incStats and ClusterSweep_incStatsRL
In this case of the last two users may be required to write their own template by editting the existing ones, see functions for details

WHAT IS THERE and INSTRUCTIONS

All scripts/function should contain a brief description and instruction at the top, a brief summary of what there is is given below.
James is happy to be contacted if you want more details
 
	CLUSTER PARAMETER SWEEPING
	Use script FullClusterSweep_singledataset.m to sweep value of dmax, and see how metrics such as numbers of clusters and composition vary.
	Also implements some assement of quality of clustering parameters as desribed in presentations to group meeting which are also avaialbe in the shared folder
	
	MICSELLANEOUS ANALYSIS FUNCTIONS
	bulkcounts - gives bulk sample counts
	COI_clustercomp - uses clusterstats file from posgen is generate cluster compositions
	clustercomp_from_indxclrpos and COI_clustercomp2, combination of these functions allows extraction of composition directly from posfiles
	clustercompvssizeplot - plot cluster compostion against a range of available size metrics. Puts clusters into size bins
	rggenerator - calcualtes radius of gyration of clusters
	edgeclusteridentifier - identifies cluster on edge of dataset via alpha hull. Also uses alpha hull to estimate dataset volume
	
	
	MORPHOLOGY AND ORIENTATION FUNCTIONS
	Combined_morpholy_DEMO.m - this script combines most of the functions below and is the "oven ready" solution, it shows how they work but more is possible.
	Fit_principle_axes - does a principle componet analysis to work out principle axes of each cluster
	Fit_princile_axes_subset - does a the same as above but only for a subset of clusters
	morphologygraph - plots "paint splatter" graph of ppt shapes
	clipCOI - allows you to determine IDS of cluster that fit a less than or more than cut off, i.e. all cluster larger than 30 or all clusters that are more oblate than prolate
	rotateaxes - rotates principle axes found by Fit_principle_axes using a rotation matrix
	stereographplot - Plots a stereographic projection of principle axes directions
	Transformaxes - Takes xy coordinates and hkl values of poles in detectorspcae, creates, xyz vectors of the poles in reconstruction space and then finds rotation matrix to tranform between reconstruction space and crystal space
	

PS Copule of abreviations I use;
RL - relabelled
COI - cluster of interest
COIIID - cluster of interest ID, i.e. the number it has in an indexed cluster pos file


