# TODO: Would be better to create a class object to contain the parameters required, this will help reduce code
#  breakage if interfaces change.
import posgenpy
from typing import List
from lxml import etree

string_vector = List[str]


def write_xml_with_relabelling(xmlFileName: str,
                               posFile: str,
                               rangeFile: str,
                               coreIons: string_vector,
                               bulkIons: string_vector,
                               relabelled_runs: int,
                               destination_folder="",
                               dclassify=0.0,
                               knn=1,
                               dmax=0.5,
                               dbulk=0.2,
                               derode=0.2,
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
                               clusterIDPosFile="clusterID.pos"):
    """
    # TODO write some documentation
    
    Function to write an XML file to perform max-sep clustering on a pos file.

    Then, runs the same analysis on relabelled data X times and saves it in a systematic way for further analysis.

    POS files for relabelled data are not saved.

    To keep it clean, save all the files in the designated_folder.
    
    xmlFileName, posFile, rangeFile, coreIons, bulkIons, relabelled_runs are required
    
    All inputs should be strings, except the core/bulk ion lists and massRandomRelabel (boolean)
    
    Set clusterstatsFile="" to not produce this file, similar for other output files
    
    See posgen manual for explanation of terms used in the XML file.
    
    Requires lxml to run. Returns the lxml.etree._ElementTree object used to write the xml file.
    
    """

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
                                         file=str(f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{clusterstatsFile}")))
        # not-clustered stats options
        if unclusterstatsFile:
            cluster.append(etree.Element("unclusterstats",
                                         file=f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{unclusterstatsFile}"))

        if sizedistFile:
            cluster.append(
                etree.Element("sizedist", file=f"{destination_folder}{xmlFileName}_random_{random_run + 1}_{sizedistFile}"))

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
               doctype="<!DOCTYPE posscript SYSTEM \"posscript.dtd\">",
               pretty_print=True)

    return tree
