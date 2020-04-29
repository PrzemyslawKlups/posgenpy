def writeClusterXML(xmlFileName, posFile, rangeFile, coreIons, bulkIons,
                          massRandomRelabel = False,
                          dclassify = "0.0", 
                          knn="1", 
                          dmax = "0.5", 
                          dbulk="0.2", 
                          derode="0.2",
                          nmin="2",
                          nmax="-1",
                          includeUnrangedPos="true",
                          includeUnrangedStats="true",
                          clusterstatsCore="true",
                          clusterstatsBulk="true", 
                          clusterstatsPercluster="true", 
                          clusterstatsFile="cluster-stats.txt",
                          unclusterstatsFile = "unclustered-stats.txt",
                          sizedistFile="sizedist.txt",
                          clusteredPosFile="cluster.pos",
                          unclusteredPosFile="unclustered.pos",
                          clusterIDPosFile="clusterID.pos"):
    
    """
    # TODO write some documentation
    
    Function to write an XML file to perform max-sep clustering on a pos file.
    
    xmlFileName, posFile, rangeFile, coreIons, bulkIons are required
    
    All inputs should be strings, except the core/bulk ion lists and massRandomRelabel (boolean)
    
    Set clusterstatsFile="" to not produce this file, similar for other output files
    
    See posgen manual for explaination of terms used in the XML file.
    
    Requires lxml to run. Returns the lxml.etree._ElementTree object used to write the xml file.
    
    """
    
    # Input error checking
    # TODO: not complete!
    if (includeUnrangedPos.lower() != "true" and includeUnrangedPos.lower() != "false"):
        raise Exception('includeUnrangedPos should be "true" or "false", the value was: {}'.format(includeUnrangedPos))
    
    from lxml import etree
    
    root = etree.Element("posscript")
    root.append( etree.Element("version") )
    # insert pos file path
    root.append( etree.Element("posload",file=posFile))
    
    # if mass randomisation is required
    if massRandomRelabel:
        root.append( etree.Element("relabel") )
    
	# TODO - this section should in a sub-function so the same
	#	code can be called as part of cluster sweep, or another
	#	function which write multiple <cluster/> operations in a 
	#	single XML file
	
    cluster = etree.SubElement(root, "cluster")
    algorithm = etree.SubElement(cluster, "algorithm", value="maxsep")
    algorithm.append( etree.Element("dclassify", value=dclassify, knn=knn) )
    algorithm.append( etree.Element("dmax", value=dmax) )
    algorithm.append( etree.Element("dbulk", value=dbulk) )
    algorithm.append( etree.Element("derode", value=derode) )
    
    # insert range file path
    cluster.append( etree.Element("range",file=rangeFile) )

    core = etree.SubElement(cluster, "core")
    typelist = etree.SubElement(core, "typelist")
    
    # List of core ions from coreList
    for coreIon in coreIons:
        typelist.append( etree.Element("atomtype", symbol=coreIon))

    bulk = etree.SubElement(cluster, "bulk")
    typelist = etree.SubElement(bulk, "typelist")
    for bulkIon in bulkIons:
        typelist.append( etree.Element("atomtype", symbol=bulkIon))

    cluster.append( etree.Element("sizeclip", nmin=nmin, nmax=nmax))

    # switch if unranged ions are required in either/both of the POS output or stats
    if includeUnrangedPos == "true" or includeUnrangedStats == "true":
            cluster.append( etree.Element("unranged", foroutput=includeUnrangedPos, forstats=includeUnrangedStats))
    
    # cluster stats options        
    if clusterstatsFile != "":
        cluster.append( etree.Element("clusterstats", 
                                      core=clusterstatsCore, 
                                      bulk=clusterstatsBulk, 
                                      percluster=clusterstatsPercluster, 
                                      file=clusterstatsFile))
    # not-clustered stats options  
    if unclusterstatsFile != "":
        cluster.append( etree.Element("unclusterstats", file=unclusterstatsFile))
    
    if sizedistFile != "":
        cluster.append( etree.Element("sizedist", file=sizedistFile))
    
    if clusteredPosFile != "":
        cluster.append( etree.Element("clustered-pos", file=clusteredPosFile, retain="true"))
        
    if unclusteredPosFile != "":
        cluster.append( etree.Element("unclustered-pos", file=unclusteredPosFile, retain="true"))

    if clusterIDPosFile != "":
        cluster.append( etree.Element("clusterid", file=clusterIDPosFile, offset="1"))

    # pre-view with print
    #print((etree.tostring(root, pretty_print=True, encoding='utf8').decode('utf8')))

    # write XML file
    tree = root.getroottree()
    tree.write(xmlFileName, 
           doctype = "<!DOCTYPE posscript SYSTEM \"posscript.dtd\">",
              pretty_print=True)
    return tree