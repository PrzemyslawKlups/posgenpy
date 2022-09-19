# TODO: Would be better to create a class object to contain the swept_parameters required, this will help reduce code
#  breakage if interfaces change.
from typing import List
from lxml import etree
import pandas
import numpy
import matplotlib.pyplot as plt
import os
# from sklearn.neighbors import KDTree
import pandas as pd

string_vector = List[str]
integer_vector = List[int]
float_vector = List[float]


def write_xml_with_relabelling(
        xmlFileName: str,
        posFile: str,
        rangeFile: str,
        coreIons: string_vector,
        bulkIons: string_vector,
        relabelled_runs: int,
        destination_folder="",
        dclassify=0.0,
        knn=1,
        dmax=0.5,
        dbulk="dmax/2",
        derode="dmax/2",
        nminV=2,
        nmaxV=-1,
        includeUnrangedPos=True,
        includeUnrangedStats=True,
        clusterstatsCore=True,
        clusterstatsBulk=True,
        clusterstatsPercluster=True,
        clusterstatsFile="cluster-stats.txt",
        unclusterstatsFile="unclustered-stats.txt",
        sizedistFile="sizedist.txt",
        clusteredPosFile="cluster.pos",
        unclusteredPosFile="unclustered.pos",
        clusterIDPosFile="clusterID.pos",
        dtd_file_location=""
):
    """Function to write an XML file to perform max-sep clustering on a pos file.
    Then, runs the same analysis on relabelled data X times and saves it in a 
    systematic way for further analysis.
    POS files for relabelled data are not saved.
    To keep it clean, save all the files in the designated_folder.
    xmlFileName, posFile, rangeFile, coreIons, bulkIons, relabelled_runs are 
    required. All inputs should be strings, except the core/bulk ion lists 
    and massRandomRelabel (boolean).
    Set clusterstatsFile="" to not produce this file, similar for other 
    output files.
    See posgen manual for explanation of terms used in the XML file.
    Requires lxml to run. Returns the lxml.etree._ElementTree object used 
    to write the xml file.

    Args:
        xmlFileName (str): name of the xml file you want to create
        posFile (str): path to the pos file
        rangeFile (str): path to the range file
        coreIons (string_vector): list of core ions (elements) you want to 
        search for
        bulkIons (string_vector): list of bulk ions you want to exclude from 
        the cluster search 
        relabelled_runs (int): how many relabelled runs you want run, the more 
        the higher accuracy of real to random cluster ratio
        destination_folder (str, optional): dir path to where you want to save 
        the data. Defaults to "".
        dclassify (float, optional): The value specifies the "classification" 
        distance, as described by Stephenson. To disable this step, set the value 
        to zero. This parameter can be used to refine which points may be 
        utilised as clustering "core" points.. Defaults to 0.0.
        knn (int, optional): This parameter can be used to refine which points 
        may be utilised as clustering "core" points. The knn attribute specifies 
        that in order to be allowed, this k-th nearest neighbouring point must 
        be found within the core radius. If not, it cannot be used as a "core" 
        clustering point Defaults to 1.
        dmax (float, optional): This value specifies the distance within which 
        points that are nominated
        to be "core" points may be found, in order to be considered part of the same 
        cluster.. Defaults to 0.5.
        dbulk (float, optional): distance from a cluster within which "bulk"
        points, can be considered to be part of that cluster. Defaults to 0.2.
        derode (float, optional): distance from non-clustered material
        that bulk points otherwise associated to clusters should be rejected from 
        clustering. Defaults to 0.2.
        nminV (int, optional): The minimum allowable size of a cluster.
        Clusters less than this size will be rejected from the analysis. Both nmin 
        and nmax are optional,
        and one or the other can be omitted. Defaults to 2.
        nmaxV (int, optional): The minimum allowable size of a cluster.
        Clusters less than this size will be rejected from the analysis. Both nmin 
        and nmax are optional,
        and one or the other can be omitted. Defaults to -1.
        includeUnrangedPos (bool, optional): if False, it will not create this 
        file. Defaults to True.
        includeUnrangedStats (bool, optional): if False, it will not create 
        this file. Defaults to True.
        clusterstatsCore (bool, optional): if False, it will not create this 
        file. Defaults to True.
        clusterstatsBulk (bool, optional): if False, it will not create this 
        file. Defaults to True.
        clusterstatsPercluster (bool, optional): if False, it will not create 
        this file. Defaults to True.
        clusterstatsFile (str, optional): name of the file. Defaults to 
        "cluster-stats.txt".
        unclusterstatsFile (str, optional): name of the file. Defaults to 
        "unclustered-stats.txt".
        sizedistFile (str, optional): name of the file. Defaults to 
        "sizedist.txt".
        clusteredPosFile (str, optional): name of the file. Defaults to 
        "cluster.pos".
        unclusteredPosFile (str, optional): name of the file. Defaults to 
        "unclustered.pos".
        clusterIDPosFile (str, optional): name of the file. Defaults to 
        "clusterID.pos".
        dtd_file_location (str, optional): location of the dtd_file in your 
        WSL system. Defaults to "".

    Returns:
        tree: xml file
    """
    # change dbulk and derode default values and link them to be dmax/2
    if dbulk == "dmax/2":
        dbulk = dmax / 2

    if derode == "dmax/2":
        derode = dmax / 2

    root = etree.Element("posscript")
    root.append(etree.Element("version"))

    # insert pos file path
    root.append(etree.Element("posload", file=posFile))

    # TODO - this section should in a sub-function so the same
    #	code can be called as part of cluster sweep, or another
    #	function which write multiple <cluster/> operations in a
    #	single XML file

    cluster = etree.SubElement(root, "cluster")
    algorithm = etree.SubElement(cluster, "algorithm", value="maxsep")
    algorithm.append(etree.Element("dclassify", value=str(dclassify), knn=str(knn)))
    algorithm.append(etree.Element("dmax", value=str(dmax)))
    algorithm.append(etree.Element("dbulk", value=str(dbulk)))
    algorithm.append(etree.Element("derode", value=str(derode)))

    # insert range file path
    cluster.append(etree.Element("range", file=rangeFile))

    core = etree.SubElement(cluster, "core")
    typelist = etree.SubElement(core, "typelist")

    # List of core ions from coreList
    for coreIon in coreIons:
        typelist.append(etree.Element("atomtype", symbol=coreIon))

    bulk = etree.SubElement(cluster, "bulk")
    typelist = etree.SubElement(bulk, "typelist")
    for bulkIon in bulkIons:
        typelist.append(etree.Element("atomtype", symbol=bulkIon))

    cluster.append(etree.Element("sizeclip", nmin=str(nminV), nmax=str(nmaxV)))

    # switch if unranged ions are required in either/both of the POS output or stats
    if includeUnrangedPos or includeUnrangedStats:
        cluster.append(etree.Element("unranged", foroutput=str(includeUnrangedPos), forstats=str(includeUnrangedStats)))

    # cluster stats options
    if clusterstatsFile:
        cluster.append(etree.Element("clusterstats",
                                     core=str(clusterstatsCore),
                                     bulk=str(clusterstatsBulk),
                                     percluster=str(clusterstatsPercluster),
                                     file=f"{destination_folder}{xmlFileName}_{clusterstatsFile}"))
    # not-clustered stats options
    if unclusterstatsFile:
        cluster.append(etree.Element("unclusterstats",
                                     file=f"{destination_folder}{xmlFileName}_{unclusterstatsFile}"))

    if sizedistFile:
        cluster.append(etree.Element("sizedist",
                                     file=f"{destination_folder}{xmlFileName}_{sizedistFile}"))

    if clusteredPosFile:
        cluster.append(etree.Element("clustered-pos",
                                     file=f"{destination_folder}{xmlFileName}_{clusteredPosFile}", retain="true"))

    if unclusteredPosFile:
        cluster.append(etree.Element("unclustered-pos",
                                     file=f"{destination_folder}{xmlFileName}_{unclusteredPosFile}", retain="true"))

    if clusterIDPosFile:
        cluster.append(etree.Element("clusterid",
                                     file=f"{destination_folder}{xmlFileName}_{clusterIDPosFile}", offset="1"))

    # iterate X times through the relabelled runs

    for random_run in range(relabelled_runs):
        # if mass randomisation is required
        root.append(etree.Element("relabel"))

        cluster = etree.SubElement(root, "cluster")
        algorithm = etree.SubElement(cluster, "algorithm", value="maxsep")
        algorithm.append(etree.Element("dclassify", value=str(dclassify), knn=str(knn)))
        algorithm.append(etree.Element("dmax", value=str(dmax)))
        algorithm.append(etree.Element("dbulk", value=str(dbulk)))
        algorithm.append(etree.Element("derode", value=str(derode)))

        # insert range file path
        cluster.append(etree.Element("range", file=rangeFile))

        core = etree.SubElement(cluster, "core")
        typelist = etree.SubElement(core, "typelist")

        # List of core ions from coreList
        for coreIon in coreIons:
            typelist.append(etree.Element("atomtype", symbol=coreIon))

        bulk = etree.SubElement(cluster, "bulk")
        typelist = etree.SubElement(bulk, "typelist")
        for bulkIon in bulkIons:
            typelist.append(etree.Element("atomtype", symbol=bulkIon))

        cluster.append(etree.Element("sizeclip", nmin=str(nminV), nmax=str(nmaxV)))

        # switch if unranged ions are required in either/both of the POS output or stats
        # if includeUnrangedPos or includeUnrangedStats:
        #     cluster.append(
        #         etree.Element("unranged", foroutput=str(includeUnrangedPos), forstats=str(includeUnrangedStats)))

        # cluster stats options
        if clusterstatsFile:
            cluster.append(etree.Element("clusterstats",
                                         core=str(clusterstatsCore),
                                         bulk=str(clusterstatsBulk),
                                         percluster=str(clusterstatsPercluster),
                                         file=str(
                                             f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{clusterstatsFile}")))
        # not-clustered stats options
        if unclusterstatsFile:
            cluster.append(etree.Element("unclusterstats",
                                         file=f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{unclusterstatsFile}"))

        if sizedistFile:
            cluster.append(
                etree.Element("sizedist",
                              file=f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{sizedistFile}"))

        # if clusteredPosFile:
        #     cluster.append(etree.Element("clustered-pos", file=clusteredPosFile, retain="true"))
        #
        # if unclusteredPosFile:
        #     cluster.append(etree.Element("unclustered-pos", file=unclusteredPosFile, retain="true"))
        #
        # if clusterIDPosFile:
        #     cluster.append(etree.Element("clusterid", file=clusterIDPosFile, offset="1"))

    # pre-view with print
    # print((etree.tostring(root, pretty_print=True, encoding='utf8').decode('utf8')))

    # write XML file
    tree = root.getroottree()
    tree.write(f"{destination_folder}{xmlFileName}.xml",
               doctype=f"<!DOCTYPE posscript SYSTEM \"{dtd_file_location}posscript.dtd\">",
               pretty_print=True)

    return tree


# FIXME: fix this xml_for_msm_with_relabelling.py:310: RuntimeWarning: invalid value encountered 
# in true_divide y_values = (real_local - random_local) / real_local

def calculate_cluster_size(
    df: pd.DataFrame,
    core_ions: list,
    bulk_ions: list,
    elements_to_exclude: list,
) -> pd.Series:
    """This function calculates cluster size based on elements provided.
    It can also exclude elements from the list such as Fe.

    Args:
        df (pd.DataFrame): dataframe with cluster stats from cluster-stats.txt, already decomposed 
        core_ions (list): list of core elements used to run this cluster analysis
        bulk_ions (list): list of bulk elements used to run this cluster analysis
        elements_to_exclude (list): list of elements to exclude (e.g. Fe)

    Returns:
        pd.Series: _description_
    """
    # combine all elements into one list
    elements_to_include = [e for e in core_ions+bulk_ions if e not in elements_to_exclude]

    # get cluster size series
    cluster_size = df.loc[:, elements_to_include].sum(axis=1)

    return cluster_size

def prepare_data_for_real_random_graphs(
    swept_parameters: float_vector,
    random_runs: int,
    xml_files: string_vector,
    n_min_values: integer_vector,
) -> dict:
    """
    This function is used to prepare random and real cluster data for other plotting 
    functions in this library.
    Use for files created with write_xml_with_relabelling.
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param swept_parameters: list with all the parameters you swept across your files
    :param random_runs: integer number of relabelled runs in a given cluster search
    :param n_min_values: list of n_min values you want to check in each swept_parameter
    :return: n_min values for each swept parameter that satisfies threshold fraction of 
    real clusters
    """

    # find number of clusters larger than n_min
    # find all real clusters for all d_max values
    all_data = {}

    for i, swept_parameter in enumerate(swept_parameters):

        # create a dictionary
        all_data[swept_parameter] = {
            "real": {},
            "random averages": {},
            "random values": {}
        }

        # real clusters
        local_size_dist_file = f"{xml_files[i].replace('.xml', '')}_sizedist.txt"
        local_df = pandas.read_csv(local_size_dist_file, sep="\t")

        # go through n_min values and save them as integers in the "real" key
        for n in n_min_values:
            local_df_filtered = local_df[(local_df["Cluster Size(core)"]) >= n]
            loc_size_core = local_df_filtered["Frequency "]
            loc_clusters = loc_size_core.sum()
            all_data[swept_parameter]["real"][n] = loc_clusters

        # random clusters
        # create an empty list to append all the random
        for n in n_min_values:
            all_data[swept_parameter]["random values"][n] = []

        # start from iterating through r to minimise pandas read csv function
        for random_run in range(random_runs):
            local_size_dist_file = f"{xml_files[i].replace('.xml', '')}_random_{random_run + 1}_sizedist.txt"
            try:
                local_df = pandas.read_csv(local_size_dist_file, sep="\t")
            except:
                continue

            # find the total number of clusters by summing up Size(core) column
            for n in n_min_values:
                local_df_filtered = local_df[(local_df["Cluster Size(core)"]) >= n]
                loc_size_core = local_df_filtered["Frequency "]
                loc_clusters = loc_size_core.sum()
                all_data[swept_parameter]["random values"][n].append(loc_clusters)

        # calculate the averages for each
        for n in n_min_values:
            all_data[swept_parameter]["random averages"][n] = numpy.array(
                all_data[swept_parameter]["random values"][n]).mean()

    return all_data


def plot_real_cluster_ratio_across_swept_param(
        swept_parameters: float_vector,
        swept_parameter_name: str,
        random_runs: int,
        xml_files: string_vector,
        n_min_values: integer_vector,
        threshold: float = 0.95,
        graph_alpha: float = 0.6,
        figsize: tuple = (9, 5),
        using_jupyter=True,
) -> None:
    """
    After sweeping through cluster search swept_parameter with posgen and generating 
    the results with randomised data, use this file to retrieve and plot the real to 
    random clusters distribution.
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param swept_parameters: list with all the parameters you swept across your files
    :param swept_parameter_name: name of the swept_parameter that's being swept, needed 
    for the graphs
    :param random_runs: integer number of relabelled runs in a given cluster search
    :param n_min_values: list of n_min values you want to check in each swept_parameter
    :param threshold: ratio of real clusters wanted in the analysis
    :return: n_min values for each swept parameter that satisfies threshold fraction of 
    real clusters
    """

    all_data = prepare_data_for_real_random_graphs(
        swept_parameters=swept_parameters,
        random_runs=random_runs,
        xml_files=xml_files,
        n_min_values=n_min_values
    )

    # plot the graph
    colors = plt.cm.viridis(numpy.linspace(0, 1, len(n_min_values)))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])
    # plt.subplot(211)
    for i, n in enumerate(n_min_values):

        # sort the data to get values for each n
        parameter_local = numpy.array([])
        real_local = numpy.array([])
        random_local = numpy.array([])

        for swept_parameter in all_data:
            parameter_local = numpy.append(parameter_local, swept_parameter)
            real_local = numpy.append(real_local, all_data[swept_parameter]["real"][n])
            random_local = numpy.append(random_local, all_data[swept_parameter]["random averages"][n])

        y_values = (real_local - random_local) / real_local

        # plot
        ax.plot(parameter_local, y_values, 'o-', color=colors[i], label="$N_{min}=%s$" % (n), alpha=graph_alpha)

    # adjust the plot
    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("$(N_{det}-N_{rand})/N_{det}$")
    ax.plot([min(swept_parameters), max(swept_parameters)], [threshold, threshold], 'r-', label=f"{round(threshold*100)}%")
    ax.legend(bbox_to_anchor=(1.05, 1))
    # ax.grid()
    # ax.xlim(left=0)
    ax.set_ylim(top=1.1)
    ax.set_title("Fraction of Real Clusters")

    # skip plt.show() if run in jupyter notebook
    if "JPY_PARENT_PID" in os.environ:
        pass
    elif using_jupyter:
        pass
    else:
        plt.show()

    # todo return exact data points from the graph
    return None

# TODO go through all plotting functions and make sure Order values are plotted as integers 

def plot_real_clusters_across_swept_param(
        swept_parameters: float_vector,
        swept_parameter_name: str,
        random_runs: int,
        xml_files: string_vector,
        n_min_values: integer_vector,
        volume: float = 0,
        graph_alpha: float = 0.7,
        figsize: tuple = (9, 5),
        using_jupyter=True, 
) -> None:
    """
    After sweeping through cluster search swept_parameter with posgen and generating the 
    results with randomised data, use this file to retrieve and plot the real to random 
    clusters distribution.
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param swept_parameters: list with all the parameters you swept across your files
    :param swept_parameter_name: name of the swept_parameter that's being swept, needed 
    for the graphs
    :param random_runs: integer number of relabelled runs in a given cluster search
    :param n_min_values: list of n_min values you want to check in each swept_parameter
    :param volume: float volume of the dataset in m^3
    :param graph_alpha: float opaqueness of the lines and markers
    :return: n_min values for each swept parameter that satisfies threshold fraction of 
    real clusters
    """

    all_data = prepare_data_for_real_random_graphs(
        swept_parameters=swept_parameters,
        random_runs=random_runs,
        xml_files=xml_files,
        n_min_values=n_min_values
    )

    # plot the same graph but with total number of clusters
    colors = plt.cm.viridis(numpy.linspace(0, 1, len(n_min_values)))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    # TODO: find the volume of the given tip to show no. density instead
    show_volume = False
    if volume != 0:
        show_volume = True

    for i, n in enumerate(n_min_values):

        # sort the data to get values for each n
        parameter_local = numpy.array([])
        real_local = numpy.array([])
        # random_local = numpy.array([])
        if show_volume:
            no_density_local = numpy.array([])

        for swept_parameter in all_data:
            parameter_local = numpy.append(parameter_local, swept_parameter)
            real_local = numpy.append(real_local, all_data[swept_parameter]["real"][n])
            # random_local = numpy.append(random_local, all_data[swept_parameter]["random averages"][n])
            if show_volume:
                no_density_local = real_local / volume
        #     y_values = (real_local - random_local)/real_local

        # plot
        if show_volume:
            ax.plot(parameter_local, no_density_local, 'o-', color=colors[i], label="$N_{min}=%s$" % (n), alpha=graph_alpha)            
        else:
            ax.plot(parameter_local, real_local, 'o-', color=colors[i], label="$N_{min}=%s$" % (n), alpha=graph_alpha)

    # adjust the plot
    ax.set_xlabel(swept_parameter_name)
    if show_volume:
        ax.set_ylabel("Number density [$m^{-3}$]")        
    else:
        ax.set_ylabel("Number of identified clusters")
    ax.legend(bbox_to_anchor=(1.3, 1))
    # ax.grid()
    # ax.set_xlim(left=0)
    # ax.set_ylim(top=1.1)
    ax.set_title("Real Clusters and %s Values" % swept_parameter_name)

    # skip plt.show() if run in jupyter notebook
    if "JPY_PARENT_PID" in os.environ:
        pass
    elif using_jupyter:
        pass
    else:
        plt.show()

    # todo return exact data points from the graph
    return None


def find_smallest_nmin(
        swept_parameters: float_vector,
        swept_parameter_name: str,
        random_runs: int,
        xml_files: string_vector,
        n_min_values: integer_vector,
        threshold: float = 0.95,
        return_real_clusters=False,
) -> integer_vector:

    """
    Reads real and random cluster searches, finds the real cluster ratio for multiple n_values
    and returns a list of smallest n_min values for each swept_parameter that satisfies given threshold
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param swept_parameters: list with all the parameters you swept across your files
    :param swept_parameter_name: name of the swept_parameter that's being swept, needed for the graphs
    :param random_runs: integer number of relabelled runs in a given cluster search
    :param n_min_values: list of n_min values you want to check in each swept_parameter
    :param threshold: fraction of real clusters wanted in the analysis
    :return: n_min values for each swept parameter that satisfies threshold fraction of real clusters
    """

    all_data = prepare_data_for_real_random_graphs(
        swept_parameters=swept_parameters,
        random_runs=random_runs,
        xml_files=xml_files,
        n_min_values=n_min_values
    )

    # find the smallest n_min in all_data:
    n_min_final = []
    real_clusters = []
    random_clusters = []

    for i, swept_parameter in enumerate(swept_parameters):
        
        n_min_final.append(0)
        real_clusters.append(all_data[swept_parameter]["real"][n_min_values[0]])
        random_clusters.append(all_data[swept_parameter]["random averages"][n_min_values[0]])

        for n in n_min_values:
            local_real = all_data[swept_parameter]["real"][n]
            # try:
            local_random = all_data[swept_parameter]["random averages"][n]
            # except:
            #     local_random = 0
            # try:
            fraction_of_real_clusters = (local_real - local_random) / local_real
            # except:
            #     fraction_of_real_clusters = 0

            if fraction_of_real_clusters > threshold:
                n_min_final[i] = int(n)
                # add local real clusters to plot them later - we want to maximise it after all
                # TODO update the code in the repo
                real_clusters[i] = int(local_real)
                random_clusters[i] = int(local_random)
                print(f"{swept_parameter_name}: {swept_parameter} | n_min: {n} | fraction: {round(fraction_of_real_clusters, 3)} | real clusters: {local_real}")
                break

    if return_real_clusters:
        return n_min_final, real_clusters, random_clusters
    else:
        return n_min_final


def prepare_data_for_composition_graphs(
        swept_parameters: float_vector,
        core_ions: string_vector,
        xml_files: string_vector,
        n_min_values_for_swept_parameters: integer_vector,
        exclude_ions=None
) -> dict:
    """
    Prepare your real cluster searches for composition plots
    :param swept_parameters:
    :param core_ions:
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param n_min_values_for_swept_parameters:
    :param exclude_ions:
    :param overlap_ions:
    :return: corrected_compositions
    """

    if exclude_ions is None:
        exclude_ions = []
    # create a new dictionary every time you run a correction to avoid correcting data more than once
    corrected_compositions = {}

    for i, swept_parameter in enumerate(swept_parameters):
        # open the cluster stats file and create pandas dataframe
        cluster_stats_file = xml_files[i].replace(".xml", "_cluster-stats.txt")
        local_df = pandas.read_csv(cluster_stats_file, sep="\t")

        # insert all the corrections in here
        # exclude (matrix) ions by deleting the column
        local_df = local_df.drop(exclude_ions, axis=1)

        # todo add this functionality later if needed
        # # lower the amount of ions for each element that has overlap with matrix ions
        # for overlap_ion in overlap_ions:
        #     local_df[overlap_ion] = round(local_df[overlap_ion] * overlap_ions[overlap_ion])
        # # transfer the overlapping ions into "original" elements
        # local_df["Ni"] = local_df["Ni"] + local_df["Fm"]
        # # delete the overlaping ions by multiplying by 0 to keep the total amount of ions constant
        # local_df["Fm"] = local_df["Fm"] * 0
        # add total ions in cluster value

        # local_df["Cluster Size"] = local_df.iloc[:, 3:-2].sum(axis=1) # no longer needed
        
        # change n_min
        local_df = local_df[local_df["Cluster Size"] >= n_min_values_for_swept_parameters[i]]

        # iterate through the core elements and import sums of the ions detected
        corrected_compositions[swept_parameter] = {}
        for ion in core_ions:
            corrected_compositions[swept_parameter][ion] = numpy.array(local_df[ion]).sum()

        # find total sum of all ions detected in the cluster search for the uncorrected composition
        corrected_compositions[swept_parameter]["total clusters"] = numpy.array(local_df.iloc[:, 3:-2]).sum()

    return corrected_compositions


def plot_cluster_composition_across_swept_param_absolute(
        swept_parameters: float_vector,
        core_ions: string_vector,
        xml_files: string_vector,
        n_min_values_for_swept_parameters: integer_vector,
        swept_parameter_name: str,
        exclude_ions=None,
        graph_alpha: float = 0.7,
        ion_colors: dict = {},
        figsize: tuple = (9, 5),
        using_jupyter=True,
) -> None:
    """
    Plot graph with absolute composition of core ions across swept parameter values.
    :param swept_parameters:
    :param core_ions:
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param n_min_values_for_swept_parameters:
    :param swept_parameter_name:
    :param exclude_ions:
    :param overlap_ions:
    :return:
    """

    # prepare the data
    corrected_compositions = prepare_data_for_composition_graphs(
        swept_parameters=swept_parameters,
        core_ions=core_ions,
        xml_files=xml_files,
        n_min_values_for_swept_parameters=n_min_values_for_swept_parameters,
        exclude_ions=exclude_ions
    )

    # plot using the data prepared
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    for ion in core_ions:
        x_values = []
        y_values = []
        for swept_parameter in swept_parameters:
            x_values.append(swept_parameter)
            try:
                y_values.append(corrected_compositions[swept_parameter][ion] /
                                corrected_compositions[swept_parameter]["total clusters"] * 100)
            except ZeroDivisionError:
                y_values.append(0)
        if ion_colors:
            ax.plot(x_values, y_values, 'o-', alpha=graph_alpha, label=ion, color=ion_colors[ion])
        else:
            ax.plot(x_values, y_values, 'o-', alpha=graph_alpha, label=ion)

    # adjust the plot
    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Cluster composition [at.%]")
    ax.legend(bbox_to_anchor=(1.2, 1))
    # ax.grid()
    # ax.set_xlim(left=0)
    # ax.set_ylim(top=1.1)
    ax.set_title(f"Absolute Cluster Composition and %s Values" % swept_parameter_name)

    # skip plt.show() if run in jupyter notebook
    if "JPY_PARENT_PID" in os.environ:
        pass
    elif using_jupyter:
        pass
    else:
        plt.show()


    return None


def plot_cluster_composition_across_swept_param_relative(
        swept_parameters: float_vector,
        swept_parameter_name: str,
        core_ions: string_vector,
        xml_files: string_vector,
        n_min_values_for_swept_parameters: integer_vector,
        exclude_ions=None,
        graph_alpha: float = 1,
        ion_colors: dict = {},
        figsize: tuple = (9, 5),
        using_jupyter=True,
) -> None:

    """
    Plot stacked columns graphs with relative composition of core ions so that sum 
    of all core ion composition = 100%
    :param swept_parameters: list of swept values
    :param swept_parameter_name: name of the parameter that is being swept
    :param core_ions: core ions you want to plot
    :param xml_files: list of full path strings with xml files for posgen cluster search
    :param n_min_values_for_swept_parameters: selected nmin values for each swept parameter
    :param exclude_ions: don't include this element in the calculations
    :return:
    """

    # prepare the data
    corrected_compositions = prepare_data_for_composition_graphs(
        swept_parameters=swept_parameters,
        core_ions=core_ions,
        xml_files=xml_files,
        n_min_values_for_swept_parameters=n_min_values_for_swept_parameters,
        exclude_ions=exclude_ions
    )

    # find the corrected value (normalised core ion composition)
    # create dictionary with place for d_max and list of core ion counts
    core_ions_sum = {}
    for swept_parameter in swept_parameters:
        core_ions_sum[swept_parameter] = []

    # iterate through the compositions data and append all ion counts into the list
    for ion in core_ions:
        for i, swept_parameter in enumerate(swept_parameters):
            core_ions_sum[swept_parameter].append(corrected_compositions[swept_parameter][ion])

    # find the total sum of all core ions for a given d_max
    for swept_parameter in core_ions_sum:
        core_ions_sum[swept_parameter] = numpy.array(core_ions_sum[swept_parameter]).sum()

    # the same graph but stacked bars
    # calculate relative composition of core ions in the clusters and plot them
    next_values = None
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.15, 0.3, 0.5, 0.6])

    # set a clear width of the columns depending on the swept parameters
    diff = numpy.zeros(len(swept_parameters)-1)
    if len(swept_parameters) > 1:
        for i in range(len(swept_parameters)-1):
            diff[i] = numpy.abs(swept_parameters[i] - swept_parameters[i+1])
        column_width = diff.mean() * 0.75
    else:
        column_width = 0.75

    for i, ion in enumerate(core_ions):
        x_values = []
        y_values = []
        for swept_parameter in swept_parameters:
            x_values.append(swept_parameter)
            y_values.append(corrected_compositions[swept_parameter][ion] / core_ions_sum[swept_parameter] * 100)

        # need to plot bar graphs on top of each other starting from 100%
        # and subtracting the values from the previous one
        # I use next value to avoid problems with labelling
        if i == 0:
            y_values_bar = numpy.ones(len(swept_parameters)) * 100

        else:
            y_values_bar = next_values

        if ion_colors:
            ax.bar(x_values, y_values_bar, alpha=graph_alpha, label=ion, width=column_width, color=ion_colors[ion])
        else:
            ax.bar(x_values, y_values_bar, alpha=graph_alpha, label=ion, width=column_width)
        next_values = y_values_bar - y_values

    # adjust the plot
    ax.set_xlabel(swept_parameter_name)
    ax.set_ylabel("Corrected relative cluster composition [at.%]")
    ax.legend(bbox_to_anchor=(1.05, 1))
    # # ax.grid()
    # ax.set_xlim(left=0)
    # ax.set_ylim(top=1.1)
    ax.set_title(f"Relative Cluster Composition and %s Values" % swept_parameter_name)

    # skip plt.show() if run in jupyter notebook
    if "JPY_PARENT_PID" in os.environ:
        pass
    elif using_jupyter:
        pass
    else:
        plt.show()

    return None

def seconds_to_hhmmss(seconds:float):
    hh = int(numpy.floor(seconds/3600))
    hh_rem = seconds % 3600
    mm = int(numpy.floor(hh_rem/60))
    ss = int(numpy.floor(seconds%60))
    time_string = f"{hh}h {mm}m {ss}s"
    return time_string

# size distribution vs order
# works correctly on decomposed cluster stats file, 
# include all ions in the txt file

constant_cluster_stats_columns = [
    "X", "Y", "Z", "Unranged", "r_gyration", "Da", "Cluster ID", "x", "y", "z", 
    "Closest Cluster ID", "Closest Cluster Distance", "convexhull_volume", "convexhull_area"
]
