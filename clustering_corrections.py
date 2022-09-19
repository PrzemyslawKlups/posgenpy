# from copyreg import remove_extension
# from selectors import DefaultSelector
import pandas as pd
import struct
# import time
import plotly.express as px
import plotly.graph_objects as go
# from pyrsistent import pdeque
from tqdm import tqdm
from scipy.spatial import ConvexHull, Delaunay
import alphashape
import matplotlib.pyplot as plt
# from descartes import PolygonPatch
import numpy as np
# import plotly.graph_objects as go
# from sklearn.neighbors import KDTree
from sklearn.neighbors import NearestNeighbors


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

# pos = read_pos("/mnt/d/RTE_21_4245/pos_files/R33_31671_F.POS")

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

# f = "/mnt/d/RTE_21_4245/RTE214245_master_range_3.RRNG"
# Ions, Ranges = read_rrng(f)


def read_apt_range_file(f):
    
    this_range_file = pd.read_csv(f)
    # example or a row below
    # Range51=71.829 71.969  Name:FeO Fe:1 O:1 Color:FF0000

    not_range_part = True
    ranges_df = pd.DataFrame(columns=['number', 'lower', 'upper', 'name', 'comp', 'colour'])

    for row in this_range_file['[Ions]']:

        # start from ranges
        if "[Ranges]" in row:
            not_range_part = False
            continue
        if not_range_part:
            continue
        elif 'Number=' in row:
            continue

        # populate the df
        split_row = row.split(sep=' ')
        split_row = [x for x in split_row if x!='']

        lower = split_row.pop(0)
        number = lower.split('=')[0].replace("Range", "")
        lower = lower.split('=')[-1]
        upper = split_row.pop(0)
        name = split_row.pop(0).replace("Name:", "")
        colour = split_row.pop(-1).replace("Color:", "")
        comp = ' '.join(split_row)

        row_to_append = {
            'number': int(number),
            'lower': float(lower), 
            'upper': float(upper), 
            'name': name, 
            'comp': comp, 
            'colour': colour
            }

        ranges_df = ranges_df.append(row_to_append, ignore_index=True)

    # change the index
    ranges_df = ranges_df.set_index('number')

    return ranges_df

# Ranges = read_apt_range_file(f)

# Range pos file

def label_ions(pos,rrngs):
    """labels ions in a .pos or .epos dataframe (anything with a 'Da' column)
    with composition and colour, based on an imported .rrng file."""

    pos['comp'] = ''
    pos['colour'] = '#FFFFFF'

    for n,r in rrngs.iterrows():

        pos.loc[(pos.Da >= r.lower) & (pos.Da <= r.upper),['comp','colour']] = [r['comp'],'#' + r['colour']]

    return pos

# label_ions(pos,Ranges)
# pos['comp'] = pos.comp.str.replace(':','')

# # Create Mass Spec

# BinWidth = 0.1
# Bins = int(pos.Da.max()/BinWidth)
# pos.Da.plot.hist(bins=Bins,
#                 alpha=1)
# plt.yscale('log')
# plt.ylim(0, 0.4e6)


def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

def not_in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)<0



def remove_edge_clusters(
    ranged_pos_df:pd.DataFrame, 
    cluster_id_pos:pd.DataFrame, 
    cluster_stats_df_all:pd.DataFrame,
    nmin:int=0,
    print_log=False
    ):

    # copy files
    print("Preparing data for edge cluster removal")
    pos = ranged_pos_df.copy()
    this_cluster_stats_df = cluster_stats_df_all.copy()
    this_cluster_id_pos = cluster_id_pos.copy()
    
    # prepare variables before the loop
    cluster_id_list = this_cluster_id_pos.loc[:, 'Da'].unique()
    edge_clusters = []
    nmin_removed = 0
    cluster_id_pos_grouped = this_cluster_id_pos.groupby('Da').mean()
    list_of_das = cluster_id_pos_grouped.sort_values(by=['x', 'y', 'z']).index
    this_cluster_stats_df = this_cluster_stats_df.sort_values(by=['X', 'Y', 'Z'])
    # area and volume for convex hull - experimental
    convex_hull_volumes = []
    convex_hull_areas = []

    if len(this_cluster_stats_df) != len(list_of_das):
        print("Function stopped. cluster-stats and clusterID.pos file have different lengths")
    this_cluster_stats_df['Da'] = list_of_das
    
    # prepare an alpha shape for pos file
    print("Calculating alphashape (might take a few seconds)")
    points = pos.iloc[:, :3].sample(frac=0.01).to_numpy()
    array_of_tuples = map(tuple, points)
    points_3d = (list(array_of_tuples))
    alpha_shape = alphashape.alphashape(points_3d, 0.1)
    cloud_vertices = pd.DataFrame(alpha_shape.vertices, columns=['x', 'y', 'z'])
    print("Alphashape calculated")  

    # iterate through every cluster and see if it's an edge cluster
    print("Iterating through all clusters")
    # for id in tqdm(cluster_id_list):
    for id in tqdm(list_of_das):
        
        # find point coordinates for each cluster
        cluster_point_cloud = this_cluster_id_pos.loc[this_cluster_id_pos['Da']==id, ['x', 'y', 'z']]

        # skip if this cluster is less than nmin 
        if len(cluster_point_cloud) < nmin:
            nmin_removed += 1
            convex_hull_volumes.append(0)
            convex_hull_areas.append(0)
            continue

        # find convex hull around each cluster
        convex_hull = ConvexHull(cluster_point_cloud)
        cluster_shell_points = cluster_point_cloud.iloc[convex_hull.vertices]
        convex_hull_volumes.append(convex_hull.volume)
        convex_hull_areas.append(convex_hull.area)

        # check if there are any points outside hull
        is_it_edge = any(not_in_hull(cluster_shell_points, cloud_vertices))
        
        # if edge cluster
        if is_it_edge:
        
            # add it to the list
            edge_clusters.append(id)

    # add area and volume columns to the main dataframe
    this_cluster_stats_df.loc[:, "convexhull_area"] = convex_hull_areas
    this_cluster_stats_df.loc[:, "convexhull_volume"] = convex_hull_volumes

    # remove edge clusters from stat files
    cluster_stats_df_no_edge = this_cluster_stats_df.loc[~this_cluster_stats_df['Da'].isin(edge_clusters)]
    cluster_stats_df_edge_only = this_cluster_stats_df.loc[this_cluster_stats_df['Da'].isin(edge_clusters)]

    # remove it from the cluster_id_pos_file 
    this_cluster_id_pos = this_cluster_id_pos.loc[~this_cluster_id_pos['Da'].isin(edge_clusters), :]

    # print how many edge clusters were removed
    all_clusters_number = len(cluster_id_list)
    clusters_removed_number = len(edge_clusters)
    current_clusters_number = len(cluster_stats_df_no_edge)
    fraction_of_clusters_removed_pct = round(clusters_removed_number/current_clusters_number*100, 1)

    # safety check
    if (current_clusters_number + clusters_removed_number) != all_clusters_number:
        print(f"Clusters don't add up!")

    if print_log:
        print(f"Clusters smaller than nmin: {nmin_removed} out of {all_clusters_number}")
        print(f"Edge clusters: {clusters_removed_number} out of {all_clusters_number-nmin_removed} ({fraction_of_clusters_removed_pct}%)") 


    # combine it into one dictionary
    data = {
        'cluster_stats_df_no_edge': cluster_stats_df_no_edge,
        'cluster_stats_df_edge_only':cluster_stats_df_edge_only, 
        'cluster_id_pos':this_cluster_id_pos, 
        'edge_clusters':edge_clusters
        }

    return data



def is_number(s):
    """Checks if a given character is a number even if its a string

    Args:
        s (string): a string

    Returns:
        bool: True if the given string can be converted to float
    """
    try:
        float(s)
        return True
    except ValueError:
        return False



def is_uppercase(s):
    """Checks if a given string is uppercase

    Args:
        s (string): string (can have more than one character)

    Returns:
        bool: True if this arg is equal to arg.upper()
    """
    try: 
        if s == s.upper():
            if is_number(s):
                print("it's a number")
                return False
            else:
                return True
        else:
            return False
    except SyntaxError:
        print('not a string')



constant_cluster_stats_columns = [
    "X", "Y", "Z", "Unranged", "r_gyration", "Da", "Cluster ID", "x", "y", "z", 
    "Closest Cluster ID", "Closest Cluster Distance", "convexhull_volume", "convexhull_area",
    "Cluster Size", "Cluster Size with Fe"
]



def decompose_cluster_stats(cluster_stats_df:pd.DataFrame):
    """Decomposes cluster stats file from POSGEN(py). Checks each column and if it has more than one ion.
    Then it adds them up to the correct columns and deletes the multi-ion columns. 
    Works for unclustered_pos file too

    Args:
        cluster_stats_df (pd.DataFrame): cluster-stats.txt file from POSGEN

    Returns:
        pd.DataFrame: copy of the cluster_stats_df with decomposed columns
    """
    # create a copy
    new_csdf = cluster_stats_df.copy()
    # find which columns need to be iterated over
    columns_for_iteration = [col for col in new_csdf.columns if col not in constant_cluster_stats_columns]

    # start iterating through each molecule
    for col in columns_for_iteration:
        
        # skip if it's only one letter
        if len(col) == 1:
            continue

        # iterate through each col name to find atoms and numbers
        col_atoms = {}
        current_atom = str(col[0])
        current_number = 1
        for letter in col[1:]: # start from the second letter to avoid first uppercase
            if is_number(letter):
                current_number = int(letter)
            elif is_uppercase(letter): # add previous atom and number to the dictionary
                col_atoms[current_atom] = current_number
                current_atom = letter
                current_number = 1
            else: # lower case letter
                current_atom += letter
            # add the final letter (in case no number or upper case at the end)
            col_atoms[current_atom] = current_number
        
        # if only one molecule and number == 1, go to the next column:
        if list(col_atoms.values()) == [1]:
            continue
        # what if there's something to decompose
        else: 
            for key in col_atoms:
                atom_number = col_atoms[key]
                if key not in new_csdf.columns:
                    new_csdf.loc[:, key] = new_csdf.loc[:, col] * atom_number
                else:
                    new_csdf.loc[:, key] = new_csdf.loc[:, key] + (new_csdf.loc[:, col] * atom_number)
            
            # change numbers in col to zero
            new_csdf.loc[:, col] = 0

    return new_csdf
    


def deconvolute_Fm_clusters(cluster_stats_file_df:pd.DataFrame) -> pd.DataFrame:
    """Deconvolutes the 29 peak labelled as Fm and adds the atoms to Ni and Fe.
    Works out max Ni expected in Fm peak based on the rest of Ni being 32% 
    of natural abundance and sets it as a limit of Fm atoms that can be transferred
    to Ni. The remaining Fm atoms are transferred to Fe. At the end, Fm column is zeroed.

    Args:
        cluster_stats_file_df (pd.DataFrame): any cluster stats file from posgen with Fe, Fm, and Ni columns

    Returns:
        pd.DataFrame: df with corrected Ni and Fe and zeroed Fm column
    """
    # check if Fm, Ni, and Fe are in there
    for col in ['Fm', 'Ni', 'Fe']:
        if col not in cluster_stats_file_df.columns:
            print(f"{col} not in the given cluster stats file")
            return False

    # rule for Fm: find the max Fm Ni atoms 
    max_Ni_in_Fm = cluster_stats_file_df.loc[:, 'Ni'] / 0.32 * 0.68
    
    # then find how much you can add to nickel
    ni_to_be_added = max_Ni_in_Fm
    ni_to_be_added[cluster_stats_file_df.loc[:, 'Fm'] < max_Ni_in_Fm] = cluster_stats_file_df.loc[:, 'Fm']
    ni_to_be_added = ni_to_be_added.apply(lambda x: int(round(x)))

    # the rest of Fm can go to Fe
    fe_to_be_added = cluster_stats_file_df.loc[:, 'Fm'] - ni_to_be_added

    # update the columns
    new_csdf = cluster_stats_file_df.copy()
    new_csdf.loc[:, 'Fe'] = new_csdf.loc[:, 'Fe'] + fe_to_be_added
    new_csdf.loc[:, 'Ni'] = new_csdf.loc[:, 'Ni'] + ni_to_be_added
    new_csdf.loc[:, 'Fm'] = 0

    # check if old and new cluster stat files have the same total number of atoms
    columns_to_iterate = [col for col in new_csdf.columns if col not in constant_cluster_stats_columns]
    if new_csdf.loc[:, columns_to_iterate].sum().sum() != cluster_stats_file_df.loc[:, columns_to_iterate].sum().sum():
        print("Total atoms in new deconvoluted file is not equal to total sum of the input cluster stats file.")

    return new_csdf


def remove_clusters_smaller_than_nmin(
    cluster_stats_df:pd.DataFrame, unclustered_stats_df:pd.DataFrame, nmin:int
):

    """Go through cluster stats file and remove clusters that contain less than nmin atoms.
    Then add these atoms back to the unclustered stats file.

    Returns:
        tuple: returns two pd.DataFrame files, first one with corrected clusters-stats
        and second one with unclustered-stats file (matrix)
    """

    # copy both files
    cluster_stats_df_new = cluster_stats_df.copy()
    unclustered_stats_df_new = unclustered_stats_df.copy()

    # remove clusters that don't satisty Nmin
    ion_columns = [col for col in cluster_stats_df_new.columns if col not in constant_cluster_stats_columns]
    
    # find clusters with sizes smaller than the threshold
    # cluster_sizes = cluster_stats_df_new.loc[:, ion_columns].sum(axis=1) #  no longer needed
    cluster_sizes = cluster_stats_df_new['Cluster Size']
    cluster_stats_to_remove = cluster_stats_df_new.loc[cluster_sizes < nmin, :]
    
    # add them to the matrix
    unclustered_stats_df_new.loc[0, ion_columns] = unclustered_stats_df_new.loc[0, ion_columns] + cluster_stats_to_remove.sum()
    # and remove from the stats file
    cluster_stats_df_nmin_removed = cluster_stats_df_new.loc[cluster_sizes >= nmin, :]

    return cluster_stats_df_nmin_removed, unclustered_stats_df_new


def get_cluster_metrics(
    reconstruction_name:str,
    reconstruction_location:str,
    sample_type:str,
    cluster_stats_df:pd.DataFrame,
    cluster_stats_df_no_edge:pd.DataFrame,
    unclustered_stats_df:pd.DataFrame,
    core_ions:list,
    bulk_ions:list,
    dmax:float,
    order:int,
    n_min:int, 
    detector_efficiency:float=0.52,
    atomic_density_per_nm3:float=85.49,
    excluded_Fe:bool=False,
    cluster_analysis_type:str="main" # main or cv01 
)->pd.DataFrame:

    """calculates cluster analysis metrics: volume fraction, number density etc. considering edge clusters.
    Requires decomposed (and deconvoluted) cluster_stat files.
    cluster_stats_df_no_edge is an output of remove_edge_clusters from the same cluster_stats_df file 

    Returns:
        pd.DataFrame: df with one row only
    """

    # copy_necessary_cols only
    columns_to_copy = [col for col in cluster_stats_df.columns if col not in constant_cluster_stats_columns]
    cluster_stats_df_new = cluster_stats_df.loc[:, columns_to_copy].copy()

    # calculate number of atoms using cluster-stats files
    ranged_atoms_in_matrix = unclustered_stats_df.sum().sum() - unclustered_stats_df.loc[0, 'Unranged']
    ranged_atoms_in_clusters = cluster_stats_df_new.loc[:, columns_to_copy].sum().sum()
    ranged_atoms_in_bulk = ranged_atoms_in_clusters + ranged_atoms_in_matrix

    # print(ranged_atoms_in_bulk, ranged_atoms_in_clusters, ranged_atoms_in_matrix)

    # clustered volume
    clustered_volume = ranged_atoms_in_clusters / (detector_efficiency * atomic_density_per_nm3)
    # volume fraction via precipitated atoms / all atoms
    volume_fraction = ranged_atoms_in_clusters / ranged_atoms_in_bulk

    # recalculate the above variables if need to exclude Fe
    if excluded_Fe:
        Fe_in_clusters = cluster_stats_df_new.loc[:, 'Fe'].sum()
        clustered_volume_no_Fe = (ranged_atoms_in_clusters - Fe_in_clusters) / (detector_efficiency * atomic_density_per_nm3)
        volume_fraction_no_Fe = (ranged_atoms_in_clusters - Fe_in_clusters) / ranged_atoms_in_bulk

    # tip volume
    tip_volume = ranged_atoms_in_bulk / (detector_efficiency * atomic_density_per_nm3)

    # Number density including half of edge clusters
    number_of_edge_clusters = len(cluster_stats_df) - len(cluster_stats_df_no_edge)
    number_of_all_clusters = len(cluster_stats_df)
    number_density = (len(cluster_stats_df) - (0.5*number_of_edge_clusters)) / tip_volume
    number_density_error = np.sqrt((len(cluster_stats_df) - (0.5*number_of_edge_clusters))) / tip_volume

    fe_status = "Included"
    if excluded_Fe:
        fe_status = "Excluded"
    
    metrics = {
        "reconstruction name": reconstruction_name,
        "sample type": sample_type,
        "cluster analysis type": cluster_analysis_type,
        "volume fraction": volume_fraction,
        "number density per nm3": number_density,
        "number density per nm3 error": number_density_error,
        "clustered volume nm3": clustered_volume,
        "tip volume nm3": tip_volume,
        "number of all clusters": number_of_all_clusters,
        "number of non-edge clusters": number_of_all_clusters - number_of_edge_clusters,
        "number of edge clusters": number_of_edge_clusters,
        "d_max [nm]": dmax,
        "Order": order,
        "N_min": n_min,
        "core ions": "+".join(core_ions),
        "bulk ions": "+".join(bulk_ions),
        "detector efficiency": detector_efficiency,
        "assumed atomic density per nm3": atomic_density_per_nm3,
        "reconstruction_location": reconstruction_location,
    }

    if excluded_Fe:
        metrics["clustered volume no Fe nm3"] = clustered_volume_no_Fe
        metrics["volume fraction no Fe"] = volume_fraction_no_Fe
        metrics["Fe excluded?"]: fe_status

    return metrics


def correct_cluster_files(
    pos_file_path:str,
    rrng_file_path:str,
    cluster_stats_file_path:str,
    unclustered_stats_file_path:str,
    clusterID_pos_file_path:str,
    xml_file:str,
    show_mass_spec:bool=False,
    print_log:bool=True,
    # detector_efficiency:float=0.52,
    # atomic_density_per_nm3:float=85.49,
    deconvolute_fm:bool=False,
    nmin:int=0,
    # include_metrics:bool=False 
):

    """Umbrella function for all the corrections in this file. Spits out data necessary for further analysis 
    including get_cluster_metrics() method. Designed for analysing MnNiSi-rich and Cu-rich clusters in RPV steels.
    Can be applied to more material systems.

    Returns:
        dict: dictionary with self-explanatory variables in keys, values are pd.DataFrame(s)
    """

    # open pos and range files
    print("Opening pos and range files")
    pos = read_pos(pos_file_path)
    _, Ranges = read_rrng(rrng_file_path)
    cluster_id_pos = read_pos(clusterID_pos_file_path)

    # check if Ranges are empty (AP Suite format) and run second rrng reader
    if len(Ranges) == 0:
        print("Range file in APSuite format. Reading it again")
        Ranges = read_apt_range_file(rrng_file_path)

    print("Labeling ions")
    label_ions(pos,Ranges)
    pos['comp'] = pos.comp.str.replace(':','')

    # Create Mass Spec
    if show_mass_spec:
        print("Showing Mass-to-charge spectrum")
        BinWidth = 0.1
        Bins = int(pos.Da.max()/BinWidth)
        pos.Da.plot.hist(bins=Bins,
                        alpha=1)
        plt.yscale('log')
        plt.show()

    # filter the ions that are ranged only
    ranged_pos = pos[pos['comp']!='']

    # open (un)cluster(ed) stats file
    print("Opening cluster stats files")
    # cluster_stats_df = pd.read_csv(cluster_stats_file_path, sep='\t')
    unclustered_stats_df = pd.read_csv(unclustered_stats_file_path, sep='\t')
    _, cluster_stats_df, _ = label_clusters(xml_file) # labeled clusters, needed for knn cluster distance

    # correct them before any further transformations
    print("Decomposing cluster stats files")
    cluster_stats_df = decompose_cluster_stats(cluster_stats_df)
    unclustered_stats_df = decompose_cluster_stats(unclustered_stats_df)
    if deconvolute_fm:
        print("Deconvoluting Fm peak in cluster stats files")
        cluster_stats_df = deconvolute_Fm_clusters(cluster_stats_df)
        unclustered_stats_df = deconvolute_Fm_clusters(unclustered_stats_df)

    print("Performing KNN Cluster Distance analysis")
    nmin_mask = cluster_stats_df['Cluster Size'] >= nmin
    knn_columns = ['Closest Cluster ID', 'Closest Cluster Distance']
    cluster_stats_df.loc[nmin_mask, knn_columns] = find_knn_cluster_distance(cluster_stats_df.loc[nmin_mask, :]).loc[:, knn_columns]

    # remove edge clusters and calculate the metrics
    print("Identifying and removing edge clusters")
    data = remove_edge_clusters(
        ranged_pos_df=ranged_pos, 
        cluster_id_pos=cluster_id_pos, 
        cluster_stats_df_all=cluster_stats_df,
        nmin=nmin,
        print_log=print_log 
    )

    # open new cluster_stats files after edge cluster removal and apply nmin filter
    print("Applying Nmin threshold to cluster stats files")
    cluster_stats_df_no_edge, unclustered_stats_df = remove_clusters_smaller_than_nmin(
        data['cluster_stats_df_no_edge'],
        unclustered_stats_df,
        nmin=nmin
    )
    cluster_stats_df_edge_only, unclustered_stats_df = remove_clusters_smaller_than_nmin(
        data['cluster_stats_df_edge_only'],
        unclustered_stats_df,
        nmin=nmin
    )
    # apply it to both from above too
    cluster_stats_df_post_nmin_prior_edge, _ = remove_clusters_smaller_than_nmin(
        cluster_stats_df, 
        unclustered_stats_df, 
        nmin
    )

    # check if the lengths of post_edge two match with the prior_edge files
    if len(cluster_stats_df_post_nmin_prior_edge) != (len(cluster_stats_df_no_edge)+len(cluster_stats_df_edge_only)):
        print("Warning! Lengths of cluster_stats file before and after edge removal don't match after applying nmin filter.")

    # create a final dictionary
    corrected_data = {
        "cluster_stats_df_post_nmin_prior_edge": cluster_stats_df_post_nmin_prior_edge,
        "cluster_stats_df_post_nmin_post_edge": cluster_stats_df_no_edge,
        "edge_cluster_stats_df": cluster_stats_df_edge_only,
        "cluster_id_pos": data['cluster_id_pos'],
        "edge_clusters": data['edge_clusters'],
        "unclustered_stats_df": unclustered_stats_df,
    }

    print("Done")

    return corrected_data

# label cluster ids across different files
# for xml_file in xml_files:

def label_clusters(xml_file):

    """Go through cluster id pos file and:
    1) label the cluster-stats file with the same cluster IDs
    2) label the cluster.pos file with the same cluster IDs
    3) change the name of the column in clusterID.pos file from Da to Cluster ID

    Returns:
        tuple: cluster_pos_df, cluster_stats_df, cluster_id_pos_df
    """

    this_cluster_stats_file_path = xml_file.replace(".xml", "_cluster-stats.txt")
    this_cluster_id_pos_file_path = xml_file.replace(".xml", "_clusterID.pos")
    this_cluster_pos_file_path = xml_file.replace(".xml", "_cluster.pos")
    
    this_cluster_stats_df = pd.read_csv(this_cluster_stats_file_path, sep='\t')
    this_cluster_id_pos_df = read_pos(this_cluster_id_pos_file_path)
    this_cluster_pos_df = read_pos(this_cluster_pos_file_path)

    # prepare variables before the loop
    this_cluster_id_pos_df = this_cluster_id_pos_df.rename(columns={'Da':'Cluster ID'})
    cluster_id_list = this_cluster_id_pos_df.loc[:, 'Cluster ID'].unique()
    # edge_clusters = []
    # nmin_removed = 0
    cluster_id_pos_grouped = this_cluster_id_pos_df.groupby('Cluster ID').mean()
    list_of_cluster_ids = cluster_id_pos_grouped.sort_values(by=['x', 'y', 'z']).index
    this_cluster_stats_df = this_cluster_stats_df.sort_values(by=['X', 'Y', 'Z'])
    

    if len(this_cluster_stats_df) != len(list_of_cluster_ids):
        print("Function stopped. cluster-stats and clusterID.pos file have different lengths")
    
    this_cluster_stats_df['Cluster ID'] = list_of_cluster_ids

    # put cluster id in cluster pos file as an extra column
    this_cluster_pos_df = this_cluster_pos_df.sort_values(by=['x', 'y', 'z'])
    this_cluster_id_pos_df = this_cluster_id_pos_df.sort_values(by=['x', 'y', 'z'])
    this_cluster_pos_df['Cluster ID'] = this_cluster_id_pos_df['Cluster ID'].apply(lambda x: int(x))

    return this_cluster_pos_df, this_cluster_stats_df, this_cluster_id_pos_df


def find_knn_cluster_distance(df:pd.DataFrame):

    # create a copy of df
    this_df = df.copy()
    
    # create a 3 column numpy array
    xyz = this_df.loc[:, ['X','Y','Z']].to_numpy()
    # calculate nearest neighbours
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(xyz)
    # save the nearest distances and indexes of the closest points
    d, indices = nbrs.kneighbors(xyz)
    
    # indices
    knn_distances = d.T[1]
    # find the closest cluster ID
    cluster_ID_list = [id for id in this_df.loc[:, 'Cluster ID'].values]
    closest_clusterID = [cluster_ID_list[i] for i in indices.T[1]]
    
    # create or overwrite the df given
    this_df.loc[:, 'Closest Cluster Distance'] = knn_distances
    this_df.loc[:, 'Closest Cluster ID'] = closest_clusterID

    return this_df



def prepare_df_for_cluster_distr(
    xml_files:list,
    swept_parameters:list,
    swept_parameter_name:str,
    include_random:bool=False,
    n_mins:list=False
)->pd.DataFrame:

    # iterate through these files and log the data
    sd_df_list = []

    for i, xml_file in enumerate(xml_files):
        
        # import the real size dist
        # use cluster_stats_files only (frequency stats are wrong!)
        # this_cs_df = pd.read_csv(xml_file.replace(".xml", "_cluster-stats.txt"), sep='\t')
        _, this_cs_df, _ = label_clusters(xml_file) # labeled clusters
        
        col_to_iterate = [col for col in this_cs_df.columns if col not in constant_cluster_stats_columns]
        this_sd_df = this_cs_df.copy()
        # this_sd_df['Cluster Size'] = this_cs_df.loc[:, col_to_iterate].sum(axis=1) # no longer needed

        # change n_min
        if n_mins:
            this_sd_df = this_sd_df[this_sd_df["Cluster Size"] >= n_mins[i]]
        
        # sort it by size
        # this_sd_df = this_sd_df.sort_values(by="Cluster Size")
        this_sd_df = this_sd_df.sort_index()

        # save it in the cluster stats df
        this_sd_df['Type'] = 'real'
        this_sd_df[swept_parameter_name] = swept_parameters[i]

        # add the nearest neighbour distances and cluster IDs
        this_sd_df = find_knn_cluster_distance(this_sd_df)
        
        # add it to the list
        sd_df_list.append(this_sd_df)

        # skip if random clusters excluded
        if not include_random:
            continue
        
        this_cs_df = pd.read_csv(xml_file.replace(".xml", "_random_1_cluster-stats.txt"), sep='\t')
        col_to_iterate = [col for col in this_cs_df.columns if col not in constant_cluster_stats_columns]
        this_sd_df = this_cs_df.copy()
        this_sd_df['Cluster Size'] = this_cs_df.loc[:, col_to_iterate].sum(axis=1)
        if n_mins:  
            this_sd_df = this_sd_df[this_sd_df["Cluster Size"] >= n_mins[i]]
        this_sd_df['Type'] = 'random'
        this_sd_df[swept_parameter_name] = swept_parameters[i]
        sd_df_list.append(this_sd_df)

    # concat all dfs
    df = pd.concat(sd_df_list)
    
    return df

# cluster distance distribution

# def prepare_df_for_cluster_distance_distr(
#     xml_files:list,
#     swept_parameters:list,
#     swept_parameter_name:str,
#     n_mins:list,
# )->pd.DataFrame:

#     """prepare cluster stats files to plot their Cluster Distance Distribution later

#     Returns:
#         pd.DataFrame: dataframe with concatenated cluster-stats files and their CDDs 
#     """

#     # add all dfs and concatenate them later to plot as a violin graph
#     knn_distributions_list = []

#     # save locations of all cluster searches
#     cluster_locations = {}

#     # for i, cluster_stats_file in enumerate(cluster_stats_files):
#     for i, xml_file in enumerate(xml_files):
        
#         # open the cluster stats file and create pandas dataframe
#         # cluster_stats_file = xml_files[i].replace(".xml", "_cluster-stats.txt")
#         # local_df = pd.read_csv(xml_file.replace(".xml", "_cluster-stats.txt"), sep="\t")
#         _, local_df, _ = label_clusters(xml_file) # labeled clusters
        
#         # add total ions in cluster value
#         col_to_iterate = [col for col in local_df.columns if col not in constant_cluster_stats_columns]
#         local_df["Cluster Size"] = local_df.loc[:, col_to_iterate].sum(axis=1)
        
#         # change n_min
#         local_df = local_df[local_df["Cluster Size"] >= n_mins[i]]
        
#         # sort it by size
#         local_df = local_df.sort_values(by="Cluster Size")

#         # create column with negative z location 
#         local_df["Z Axis"] = local_df["Z"]

#         # NN distr
#         X = local_df[['X', 'Y', 'Z']].to_numpy()
#         if len(X) == 0:
#             continue
#         kdt = KDTree(X)
#         d, indices = kdt.query(X, k=2, return_distance=True)

#         # plot the knn connections
#         df = pd.DataFrame(X, columns=['x', 'y', 'z'])
#         indices.sort()
#         df['group'] = [float(f"{row[0]}.{row[1]}") for row in indices]

#         # save distributions of these lengths
#         this_df = pd.DataFrame()
#         this_df.loc[:, 'knn_distance'] = pd.Series(d.transpose()[1])
#         this_df[swept_parameter_name] = swept_parameters[i]
#         knn_distributions_list.append(this_df)
        
#         # save cluster locations
#         key_name = str(swept_parameters[i])
#         cluster_locations[key_name] = X

#     # combine all dfs
#     knn_df = pd.concat(knn_distributions_list)

#     return knn_df


def show_clusters_in_pos(
    clusters_to_plot:list, 
    xml_files:list,
    swept_param_values:list,
    chosen_param_value:float, 
    range_file_path:str
):

    # find selected_xml_file
    selected_xml_file = [xml_files[i] for i, x in enumerate(swept_param_values) if x == chosen_param_value]
    if len(selected_xml_file) != 1:
        print("Found 0 or 2< files. Can't continue")
        return None
    selected_xml_file = selected_xml_file[0]

    # cluster_id_pos_df = cc.read_pos(cluster_id_pos_file_path)
    cluster_id_pos_file_path = selected_xml_file.replace(".xml", "_clusterID.pos")
    cluster_id_pos_df = read_pos(cluster_id_pos_file_path)

    # find it in the cluster id and cluster data
    selected_cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df['Da'].isin(clusters_to_plot), :]

    # open pos and range files
    # print("Opening pos and range files")
    # pos = cc.read_pos(pos_file_path)
    _, Ranges = read_rrng(range_file_path)
    # cluster_id_pos = cc.read_pos(clusterID_pos_file_path)

    # check if Ranges are empty (AP Suite format) and run second rrng reader
    if len(Ranges) == 0:
        print("Range file in APSuite format. Reading it again")
        Ranges = read_apt_range_file(range_file_path)

    # label cluster pos file 
    cluster_pos = read_pos(selected_xml_file.replace(".xml", "_cluster.pos"))
    label_ions(cluster_pos, Ranges)
    cluster_pos['comp'] = cluster_pos['comp'].str.replace(':','')

    # set up the limits for the graph
    xmin = selected_cluster_id_pos_df.loc[:, 'x'].min()
    xmax = selected_cluster_id_pos_df.loc[:, 'x'].max()
    ymin = selected_cluster_id_pos_df.loc[:, 'y'].min()
    ymax = selected_cluster_id_pos_df.loc[:, 'y'].max()
    zmin = selected_cluster_id_pos_df.loc[:, 'z'].min()
    zmax = selected_cluster_id_pos_df.loc[:, 'z'].max()

    # filter the pos data
    mask1 = (cluster_pos['x'] >= xmin) & (cluster_pos['x'] <= xmax)
    mask2 = (cluster_pos['y'] >= ymin) & (cluster_pos['y'] <= ymax)
    mask3 = (cluster_pos['z'] >= zmin) & (cluster_pos['z'] <= zmax)
    all_masks = mask1 & mask2 & mask3
    pos_for_searching = cluster_pos[all_masks]

    # clean up the names
    new_names = pos_for_searching.loc[:, 'comp'].apply(lambda x: x.replace("1", "").replace("Name", ""))
    pos_for_searching.loc[:, 'comp'] = new_names
    # new_names = [name.replace("1", "").replace("Name", "") for name in comp_names]

    # create the ivas color map
    pos_for_searching_colors = pos_for_searching.loc[:, ['comp', 'colour']].drop_duplicates()
    ivas_color_map = {}
    for row in pos_for_searching_colors.to_numpy():
        ivas_color_map[row[0]] = row[1]

    # plot it in 3d
    fig = px.scatter_3d(pos_for_searching, x='x', y='y', z='z', color='comp', color_discrete_map=ivas_color_map)
    fig.update_traces(marker_size = 2)
    fig.update_layout(title=f'Cluster ID(s): {clusters_to_plot}')
    fig.show()

    return pos_for_searching



def compare_selected_clusters_across_param(
    xml_files:list,
    swept_param_values:list,
    swept_param_name:str,
    chosen_param_value:float,
    selected_clusters:list
)->go.Figure:

    
    # find selected_xml_file
    selected_xml_file = [xml_files[i] for i, x in enumerate(swept_param_values) if x == chosen_param_value]
    if len(selected_xml_file) != 1:
        print("Found 0 or 2< files. Can't continue")
        return None
    selected_xml_file = selected_xml_file[0]
    
    # find xyz range of selected clusters 
    main_pos_file_df = read_pos(selected_xml_file.replace(".xml", "_clusterID.pos"))
    main_pos_file_df = main_pos_file_df.loc[main_pos_file_df['Da'].isin(selected_clusters)]

    xmin = main_pos_file_df.x.min()
    xmax = main_pos_file_df.x.max()
    ymin = main_pos_file_df.y.min()
    ymax = main_pos_file_df.y.max()
    zmin = main_pos_file_df.z.min()
    zmax = main_pos_file_df.z.max()

    # go through each cluster id pos file 
    dfs_to_concat = []
    for i, xml_file in enumerate(xml_files):
        # open each file
        cluster_id_pos_path = xml_file.replace(".xml", "_clusterID.pos")
        cluster_id_pos_df = read_pos(cluster_id_pos_path)

        # filter pos file using xyz ranges
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["x"] >= xmin, :]
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["x"] <= xmax, :]
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["y"] >= ymin, :]
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["y"] <= ymax, :]
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["z"] >= zmin, :]
        cluster_id_pos_df = cluster_id_pos_df.loc[cluster_id_pos_df["z"] <= zmax, :]

        # assign the values for nicer plotting
        cluster_id_pos_df[swept_param_name] = str(swept_param_values[i])
        cluster_id_pos_df['Cluster ID'] = cluster_id_pos_df['Da'].apply(lambda x: str(round(x)))
        dfs_to_concat.append(cluster_id_pos_df)

    # concat all of them
    df_to_plot = pd.concat(dfs_to_concat)

    # plot 
    fig = px.scatter_3d(
        df_to_plot,
        x='x', 
        y='y',
        z='z',
        color=swept_param_name,
        symbol='Cluster ID'
    )

    fig.update_traces(marker_size=3, opacity=0.8)

    return fig

# read IVAS colours from range file
def get_element_colours_from_range_file(range_file_path:str)->dict:

    # read range file
    _, ranges = read_rrng(range_file_path)

    # massage the data
    colours_df = ranges.loc[:, ['comp', 'colour']]
    colours_df.loc[:, 'comp'] = colours_df.loc[:, 'comp'].apply(lambda x: x.replace(":1", ""))
    colours_df.loc[:, 'comp'] = colours_df.loc[:, 'comp'].apply(lambda x: x.replace("Name:", ""))
    colours_df.loc[:, 'comp'] = colours_df.loc[:, 'comp'].apply(lambda x: x.replace(":", ""))
    colours_df.loc[:, 'colour'] = colours_df.loc[:, 'colour'].apply(lambda x: "#"+x)
    colours_df = colours_df.drop_duplicates()
    colours_df = colours_df.reset_index().set_index('comp').drop(['number'], axis=1)

    # return a dictionary with elements as keys and hashed colours as values
    ivas_colour_dict = colours_df.to_dict()['colour']

    return ivas_colour_dict
