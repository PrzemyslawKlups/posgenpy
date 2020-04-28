# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 15:51:20 2019

@author: Ben
"""

# Open pos file - should be in shared folder
import pandas as pd
import struct
start = time.time()
def read_pos(f):

    """ Loads an APT .pos file as a pandas dataframe.
    Columns:
        x: Reconstructed x position
        y: Reconstructed y position
        z: Reconstructed z position
        Da: mass/charge ratio of ion"""
    # read in the data
    n = len(open(f,"rb").read())/4
    d = struct.unpack('>'+'f'*int(n),open(f,"rb").read(4*int(n)))
                    # '>' denotes 'big-endian' byte order
    # unpack data
    pos = pd.DataFrame({'x': d[0::4],
                        'y': d[1::4],
                        'z': d[2::4],
                        'Da': d[3::4]})
    return pos

pos = read_pos("C:/Users/Ben/Documents/GitHub/APT-Data-Analysis/Automated Edge Cluster Detection/Test Files/006_full.pos")

# Import range file
def read_rrng(f):
    """Loads a .rrng file produced by IVAS. Returns two dataframes of 'ions'
    and 'ranges'."""
    import re

    rf = open(f,'r').readlines()

    patterns = re.compile(r'Ion([0-9]+)=([A-Za-z0-9]+).*|Range([0-9]+)=(\d+.\d+) +(\d+.\d+) +Vol:(\d+.\d+) +([A-Za-z:0-9 ]+) +Color:([A-Z0-9]{6})')

    ions = []
    rrngs = []
    for line in rf:
        m = patterns.search(line)
        if m:
            if m.groups()[0] is not None:
                ions.append(m.groups()[:2])
            else:
                rrngs.append(m.groups()[2:])

    ions = pd.DataFrame(ions, columns=['number','name'])
    ions.set_index('number',inplace=True)
    rrngs = pd.DataFrame(rrngs, columns=['number','lower','upper','vol','comp','colour'])
    rrngs.set_index('number',inplace=True)


    rrngs[['lower','upper','vol']] = rrngs[['lower','upper','vol']].astype(float)
    rrngs[['comp','colour']] = rrngs[['comp','colour']].astype(str)

    return ions,rrngs

Ions, Ranges = read_rrng("C:/Users/Ben/Documents/Post Doc/Code/Python/006_full.cluster.rrng")

# Range pos file

def label_ions(pos,rrngs):
    """labels ions in a .pos or .epos dataframe (anything with a 'Da' column)
    with composition and colour, based on an imported .rrng file."""

    pos['comp'] = ''
    pos['colour'] = '#FFFFFF'

    for n,r in rrngs.iterrows():

        pos.loc[(pos.Da >= r.lower) & (pos.Da <= r.upper),['comp','colour']] = [r['comp'],'#' + r['colour']]

    return pos

label_ions(pos,Ranges)
pos['comp'] = pos.comp.str.replace(':','')

# Create Mass Spec

#BinWidth = 0.1
#Bins = int(pos.Da.max()/BinWidth)
#pos.Da.plot.hist(bins=Bins,
#                 alpha=1)

# Create cluster search
# Maximum separation/Core-linkage
CoreIons = {'Cu1','Ni1'}
Dmax = 1
Nmin = 40
Order = 1
L = 0.5
E = 0.5

# Filter DF to contain just core ions
CorePos = pos.loc[pos['comp'].isin(CoreIons)]

import numpy as np
import scipy.spatial as spatial

# If Core ions within Dmax is > Order, keep ions
# Want array - index of location & number within Dmax
CoreIonLocations = np.array(CorePos[['x','y','z']])
point_tree = spatial.cKDTree(CoreIonLocations)

# Link core ions that are connected
CoreIonsWithinDmax = []
for SoluteIon in range(len(CoreIonLocations)):
    CoreIonsWithinDmax.extend([(point_tree.query_ball_point(CoreIonLocations[SoluteIon], Dmax))])

# Select only ions that exceed order requirement 
# (>= Order number of core ions within Dmax of ion +1 needed to account for ion itself)
CoreIonsExceedingOrder = [Ion for Ion in CoreIonsWithinDmax if len(Ion) >= Order + 1]

# L Step prior to E step
# Array for core ions exceeding order - their x, y, z co-ordinates
# CoreIonsExceedingOrderLocations = list(set([item for sublist in CoreIonsExceedingOrder for item in sublist]))
ClusteredIonsLocations = CoreIonLocations[[list(set([item for sublist in CoreIonsExceedingOrder for item in sublist]))]]

# Create cKDTree for all ions in dataset
AllIonpoint_tree = spatial.cKDTree(np.array(pos[['x','y','z']]))

# Get ion locations of ions within L of core ions
IonsWithinLOfCoreIons = []
for ArrayLocation in range(len(ClusteredIonsLocations)):
    IonsWithinLOfCoreIons.extend([(AllIonpoint_tree.query_ball_point(ClusteredIonsLocations[ArrayLocation], L))])

ListOfIonsInClusterAfterL = list(set([item for sublist in IonsWithinLOfCoreIons for item in sublist]))
# Create list of ions in matrix after L step

# Get location atoms within E of non-clustered atoms
# Create dataframe of ions within E of non-clustered atoms
Matrixpos = pos.loc[list([set(range(len(pos))) - set(ListOfIonsInClusterAfterL)])[0]]
#
MatrixIonLocations = np.array(Matrixpos[['x','y','z']])
Matrixpoint_tree = spatial.cKDTree(np.array(Matrixpos[['x','y','z']]))
# Find list of ions within E of matrix ions
IonsWithinEOfMatrixIons = []
for ArrayLocation in range(len(MatrixIonLocations)):
    IonsWithinEOfMatrixIons.extend([(AllIonpoint_tree.query_ball_point(MatrixIonLocations[ArrayLocation], E))])

ListOfIonsWithinEOfMatrixIons = list(set([item for sublist in IonsWithinEOfMatrixIons for item in sublist]))

NonEIonLocations = np.array(pos.drop(pos.index[ListOfIonsWithinEOfMatrixIons])[['x','y','z']])
NonEIonpoint_tree = spatial.cKDTree(NonEIonLocations)

# Get ion locations of ions within L of core ions (excluding E removed atoms)
IonsWithinLOfCoreIons = []
for ArrayLocation in range(len(ClusteredIonsLocations)):
    IonsWithinLOfCoreIons.extend([(NonEIonpoint_tree.query_ball_point(ClusteredIonsLocations[ArrayLocation], L))])

ListOfIonsInClusterAfterL = list(set([item for sublist in IonsWithinLOfCoreIons for item in sublist]))

import networkx 
# from networkx.algorithms.components.connected import connected_components

def to_graph(m):
    G = networkx.Graph()
    for part in m:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(m):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(m)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

G = to_graph(IonsWithinLOfCoreIons)
Clusters = dict(list(enumerate(list(networkx.connected_components(G)))))
Clusters = { key:value for (key,value) in Clusters.items() if len(value) > Nmin}
end = time.time()
print(end - start)

# Create sweep tool to do Dmax, Nmin, Order
# Randomise ion identity and compare - subtract and calc. CIs?