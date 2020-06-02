function [] = ClusterSweep_incStatsRL(posfile,rangefile,ClusteredIons,BulkIons,dinitial,dstep,dfinal,inclusion,outputfolder)
%Read description of ClusterSweep_incStats, it's exactly the same but
%doesn't do relabel step.

%NB clusterstats will be quite affected by Nmin which isn't currently set,
%may add this option in future.
%ClusterSweep_incStats same as ClusterSweep except it also outputs all the 
%clusterstats file from posgen as well - takes pos and range file and sweeps %
%dmax using input limits, for a relabbeld dataset (currently relabelling at every value of dmax)
%Inclusion can be done, ouput is size distribution
%files + clusterstats files

%posfile = text path to desired pos file
%rangefile = text path to desired range file
%ClusteredIons = ionic species that are to be included in clustering step,
%NB current maximum is ____
%BulkIons = ionic species that are to be included in inclusion step, does
%not have to be all bulk ions%
%dinitial, dstep and dfinal = numerical values at which to start and end
%sweep inclusive and the step size
%inclusion test 'yes' if inclusion is desired. Otherwise just write random
    %number, but it doesn't like 'no'. Some input is required
%ouputfolder = text path to desired folder and name of new folder to be
    %made

    
%select correct xml input file to open, with right number of clustered and bulk ions be assigned%%could work out how to write an xml file directly, which would be nicer%
%Currently user needs to make a new xml file named in the same way if they
%required a number combination of clustered and bulk ions that doesn't%
%file name format is MatlabClusterSweep_icStats_xCI_yBI.xml  where x is
%number of cluster ions and y is no of bulk ions of interest

x = size(ClusteredIons,2);
y =size(BulkIons,2);
status = mkdir(outputfolder);
if exist(['MatlabClusterSweep_incStats_RL_' num2str(x) 'CI_' num2str(y) 'BI.xml']) == 2 %does correct input xml file exist?%
    disp('Input xml file exists')
else
    disp('Input xml file not found, please create file of correct format, see ClusterSweep_incStatsRL.m')
    return
end
DOMnode = xmlread(['MatlabClusterSweep_incStats_RL_' num2str(x) 'CI_' num2str(y) 'BI.xml']);%loads xml file which is already in correct input format%
% theStruct = parseXML('C:\Users\James\Documents\Posgen\170710-posgen\clusterSweeptest.xml')
% gives structure of xml file if needed%
for i = dinitial:dstep:dfinal %parameter sweep in dmax%
    dmax = num2str(i);
    sizedistfile = [outputfolder '\sizedist-RL-dmax-' num2str(i) '.txt'];
    clusterstatsfile = [outputfolder '\clusterstats-RL-dmax-' num2str(i) '.txt'];
    if inclusion == 'yes'
        dbulk = num2str(i/2);
    else
       dbulk ='0.0'; 
    end
%next lines edit xml file to the correct values%
DOMnode.getElementsByTagName('posload').item(0).setAttribute('file',posfile); 
DOMnode.getElementsByTagName('range').item(0).setAttribute('file',rangefile);
DOMnode.getElementsByTagName('dmax').item(0).setAttribute('value',dmax);
DOMnode.getElementsByTagName('dbulk').item(0).setAttribute('value',dbulk);
DOMnode.getElementsByTagName('derode').item(0).setAttribute('value',dbulk);%commemt out if not using%
DOMnode.getElementsByTagName('sizedist').item(0).setAttribute('file',sizedistfile);
DOMnode.getElementsByTagName('clusterstats').item(0).setAttribute('file',clusterstatsfile);

%edit names of elements%
for i = 0:x-1
DOMnode.getElementsByTagName('atomtype').item(i).setAttribute('symbol',ClusteredIons(i+1));
end
for j = 0:y-1
    DOMnode.getElementsByTagName('atomtype').item(x+j).setAttribute('symbol',BulkIons(j+1));
end

xmlwrite('posgeninput_ClusterSweepRL.xml',DOMnode);%write final xml fed to posgen%
!posgen.exe & ;
status = dos('posgen posgeninput_ClusterSweepRL.xml'); %excutes posgen%
end
end
