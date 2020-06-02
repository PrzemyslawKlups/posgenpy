%THIS SCRIPT DOES NOT SPECIFY THE CLUSTERING AND BULK ELEMENTS THIS MUST BE
%EDITTED IN THE XML MANUALLY%%posgeninputsinglecluster.xml%
function Cluster_makepos(posfile,rangefile,samplename,dmax,dbulk,derode,Nmin,outputfolder)
%Inputs%
% posfile - path to a posfile
% rangefile - path to a range file
% samplename - sample code with which to save files eg 'R5083_9999_roi'
% dmax 
% dbulk
% derode 
% Nmin
% outputfolder - folder into which the output files will be saved

%convert parameter to text%
dmax = num2str(dmax);
dbulk = num2str(dbulk);
derode = num2str(derode);
Nmin = num2str(Nmin);

disp('WARNING, THIS SCRIPT DOES NOT SPECIFY THE CLUSTERING AND BULK ELEMENTS THIS MUST BE EDITTED IN THE XML MANUALLY posgeninputsinglecluster.xml')

%Make output file names
sizedistfile = [outputfolder '\' samplename '-sizedist-dmax-' num2str(dmax) '-Nmin-' num2str(Nmin) '.txt'];
clusterstatsfile = [outputfolder '\' samplename '-clusterstats-dmax-' num2str(dmax) '-Nmin-' num2str(Nmin) '.txt'];
clusteredposfile = [outputfolder '\' samplename '-dmax-' num2str(dmax) '-Nmin-' num2str(Nmin) '.cluster.pos'];
indexclusterposfile = [outputfolder '\' samplename '-dmax-' num2str(dmax) '-Nmin-' num2str(Nmin) '.cluster.indexed.pos'];
%this could be a little more tricky if I want to keeep exactly the same
%naming conventions as IVAS%

DOMnode = xmlread(['posgeninputsinglecluster.xml']);%loads xml file which is already in correct input format%
%theStruct = parseXML('C:\Users\James\Documents\Matlab\posgen\posgeninputsinglecluster.xml');
% gives structure of xml file if needed%

%next lines edit xml file to the correct values%
DOMnode.getElementsByTagName('posload').item(0).setAttribute('file',posfile); 
DOMnode.getElementsByTagName('range').item(0).setAttribute('file',rangefile);
DOMnode.getElementsByTagName('dmax').item(0).setAttribute('value',dmax);
DOMnode.getElementsByTagName('dbulk').item(0).setAttribute('value',dbulk);
DOMnode.getElementsByTagName('derode').item(0).setAttribute('value',derode);
DOMnode.getElementsByTagName('sizedist').item(0).setAttribute('file',sizedistfile);
DOMnode.getElementsByTagName('clusterstats').item(0).setAttribute('file',clusterstatsfile);
DOMnode.getElementsByTagName('sizeclip').item(0).setAttribute('nmin',Nmin);
DOMnode.getElementsByTagName('clustered-pos').item(0).setAttribute('file',clusteredposfile);
%DOMnode.getElementsByTagName('unclustered-pos').item(0).setAttribute('file',unclusteredposfile);
DOMnode.getElementsByTagName('clusterid').item(0).setAttribute('file',indexclusterposfile);

xmlwrite('posgeninputsinglecluster.xml',DOMnode);%write final xml fed to posgen%
!posgen.exe & ;
status = dos('posgen posgeninputsinglecluster.xml'); %excutes posgen%
end