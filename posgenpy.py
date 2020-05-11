# TODO: Would be better to create a class object to contain the parameters required,
# this will help reduce code breakage if interfaces change.


def write_cluster_xml(xml_file, pos_file, range_file, core_ions, bulk_ions, *,
                      mass_random_relabel=False,
                      dclassify:float=0.0, 
                      knn: int =1,
                      dmax: float=0.5 ,
                      dbulk: float=0.2 ,
                      derode: float=0.2 ,
                      nmin: int=2 ,
                      nmax: int=-1 ,
                      include_unranged_pos=True,
                      include_unranged_stats=True,
                      clusterstats_core=True,
                      clusterstats_bulk=True,
                      clusterstats_percluster=True,
                      clusterstats_file="cluster-stats.txt",
                      unclusterstats_file="unclustered-stats.txt",
                      sizedist_file="sizedist.txt",
                      clustered_pos_file="cluster.pos",
                      unclustered_pos_file="unclustered.pos",
                      clusterid_pos_file="clusterID.pos"):
    """
    # TODO write some documentation
    
    Function to write an XML file to perform max-sep clustering on a pos file.
    
    xml_file, pos_file, range_file, core_ions, bulk_ions are required
    
    All inputs should be strings, except the core/bulk ion lists and mass_random_relabel (boolean)
    
    Set clusterstats_file="" to not produce this file, similar for other output files
    
    See posgen manual for explanation of terms used in the XML file.
    
    Requires lxml to run. Returns the lxml.etree._ElementTree object used to write the xml file.
    
    """

    from lxml import etree

    root = etree.Element("posscript")
    root.append(etree.Element("version"))
    # insert pos file path
    root.append(etree.Element("posload", file=pos_file))

    # if mass randomisation is required
    if mass_random_relabel:
        root.append(etree.Element("relabel"))

    # TODO - this section should in a sub-function so the same
    #   code can be called as part of cluster sweep, or another
    #   function which write multiple <cluster/> operations in a
    #   single XML file

    cluster = etree.SubElement(root, "cluster")
    algorithm = etree.SubElement(cluster, "algorithm", value="maxsep")
    algorithm.append(etree.Element("dclassify", value=str(dclassify), knn=str(knn)))
    algorithm.append(etree.Element("dmax", value=str(dmax)))
    algorithm.append(etree.Element("dbulk", value=str(dbulk)))
    algorithm.append(etree.Element("derode", value=str(derode)))

    # insert range file path
    cluster.append(etree.Element("range", file=range_file))

    core = etree.SubElement(cluster, "core")
    typelist = etree.SubElement(core, "typelist")

    # List of core ions from coreList
    for coreIon in core_ions:
        typelist.append(etree.Element("atomtype", symbol=coreIon))

    bulk = etree.SubElement(cluster, "bulk")
    typelist = etree.SubElement(bulk, "typelist")
    for bulkIon in bulk_ions:
        typelist.append(etree.Element("atomtype", symbol=bulkIon))

    cluster.append(etree.Element("sizeclip", nmin=str(nmin), nmax=str(nmax)))

    # switch if unranged ions are required in either/both of the POS output or stats
    if include_unranged_pos or include_unranged_stats:
        cluster.append(etree.Element("unranged",
                                     foroutput=str(include_unranged_pos),
                                     forstats=str(include_unranged_stats)))

    # cluster stats options        
    if clusterstats_file:
        cluster.append(etree.Element("clusterstats",
                                     core=str(clusterstats_core),
                                     bulk=str(clusterstats_bulk),
                                     percluster=str(clusterstats_percluster),
                                     file=str(clusterstats_file)))
    # not-clustered stats options  
    if unclusterstats_file:
        cluster.append(etree.Element("unclusterstats", file=unclusterstats_file))

    if sizedist_file:
        cluster.append(etree.Element("sizedist", file=sizedist_file))

    if clustered_pos_file:
        cluster.append(etree.Element("clustered-pos", file=clustered_pos_file, retain="true"))

    if unclustered_pos_file:
        cluster.append(etree.Element("unclustered-pos", file=unclustered_pos_file, retain="true"))

    if clusterid_pos_file:
        cluster.append(etree.Element("clusterid", file=clusterid_pos_file, offset="1"))

    # pre-view with print
    # print((etree.tostring(root, pretty_print=True, encoding='utf8').decode('utf8')))

    # write XML file
    tree = root.getroottree()
    tree.write(xml_file,
               doctype="<!DOCTYPE posscript SYSTEM \"posscript.dtd\">",
               pretty_print=True)
    return tree
