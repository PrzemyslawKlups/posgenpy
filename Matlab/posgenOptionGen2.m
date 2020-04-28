function [coreList, bulkList]=posgenOptionGen2(xmlName,posFile,rrngFile,dmax,knn,dbulk,derode,nmin,core,statFile,relabel,clustered,unclustered,clusterIndex,unranged)
% creates an options.xml (xmlName) file for use with posgen
% function posgenOptionGen(xmlName,posFile,rrngFile,dmax,knn,dbulk,derode,...
% nmin,core,statFile,FLAGS:relabel,clustered,unclustered,clusterIndex)
% core = {'Y' 'YO' 'TiO' etc} possible ions produced by rrngFile to use as
% core ions in the cluster search.
% TODO: posscriptDir input (currently hard coded ln 33) May 2016
% ONLY works with RRNG files, see rrngReader
% A London April 2013
% % *example values*
% all strings:
% xmlName = 'options.xml';
% posFile = 'POSFILE';
% rrngFile = 'r14_15065_10.rrng';
% If one of dmax, (knn), dbulk, (derode or nmin) is a vector then the
%   clustersweep command will be called in posgen and a folder of data
%   produced
% dmax = 'DMAX'; % maximum separation distance (in nm)
% knn = '1' % minimum core ions withing dmax to count as clustered, >= 1
% dbulk = 'DBULK'; % aka envelope distance (in nm)
% derode = '0'; % erosion distance
% nmin = 'NMIN'; % '12'
% core = {'Y' 'YO' 'Ti' 'TiO' 'O' 'O2'}; % potential core ions
% Other/bulk/matrix ions are identified by taking a difference of the core
%   ions from the ions defined in the range file supplied.
% statFile = 'clusterStats.txt';
% % This assumes the range file has elements in  the correct order eg Fe:1
% % O:1 for FeO not OFe.
% Strings for outpout pos-files:
% clustered = 0; % 0 for none, string for pos output
% unclustered = 0; % 0 for none, string for pos output
% clusterIndex is a pos file with masses relabelled as cluster IDs
% clusterIndex = 'cluster.indexed'; % 0 for none string for pos output
% all logic
% relabel = 0; % randomising flag
posscriptDir = 'C:\cygwin64\home\andy\posgen\';
% get all ions from rrng file using ionList
[element_num, range_num, elements, ranges] = rangeReader (rrngFile);
ions = ionList(ranges,elements);
% check core ions against list and produce coreList
coreList = ions(ismember(ions,core));
coreList = flipdim(coreList,2); % 'sorting'
% all remaining ions go in bulkList
bulkList = ions(logical(1-ismember(ions,core)));
bulkList = sort(bulkList); % sorting

sweep = 'none';
if isnumeric(dmax) && length(dmax)>1
    sweep = 'dmax';
    sweepVar = dmax;
    dmax = num2str(dmax(1));
elseif isnumeric(dbulk) && length(dbulk)>1
    sweep = 'dbulk';
    sweepVar = dbulk;
    dbulk = num2str(dbulk(1));
elseif isnumeric(derode) && length(derode)>1
    sweep = 'derode';
    sweepVar = derode;
    derode = num2str(derode(1));
end
if ~strcmp(sweep,'none');
    sweepStart = sweepVar(1);
    sweepEnd = sweepVar(end)+1E-6; % fix
    sweepStep = sweepVar(2)-sweepVar(1);
    disp(['Using sweep step of ' num2str(sweepStep) ' for ' sweep]);
end

if ~exist('unranged','var')
    unranged = false;
end

if str2double(knn)>1 % if knn > 1 then use knn order parameter
    useKnn = '1';
else
    useKnn = '0'; % else disable this clustering step
end

% write options.xml
fid = fopen(xmlName,'w');
fprintf(fid,'%s\n',strcat('<!DOCTYPE posscript SYSTEM "',regexprep(posscriptDir,'\\','/'),'posscript.dtd">'));
fprintf(fid,'%s\n','<posscript>');
fprintf(fid,'\t%s\n','<version value="0.0.1"/>');
fprintf(fid,'\t%s\n',strcat('<posload file="',posFile,'"/>'));
if(relabel)
    fprintf(fid,'\t%s\n','<relabel/>');
end
% Optional sweep tag
if ~strcmp(sweep,'none')
     fprintf(fid,'\t<clustersweep parameter="%s" start="%f" end="%f" step="%f">\n',sweep,sweepStart,sweepEnd,sweepStep);
end
fprintf(fid,'\t%s\n','<cluster>');
fprintf(fid,'\t\t%s\n','<algorithm value="maxsep">');
fprintf(fid,'\t\t\t%s\n',strcat('<dclassify value="',useKnn,'" knn="',knn,'"/> <!-- Coring (pre-cluster) dist; zero to disable this step-->'));
fprintf(fid,'\t\t\t%s\n',strcat('<dmax value="',dmax,'"/> <!-- Max sep dist -->'));
fprintf(fid,'\t\t\t%s\n',strcat('<dbulk value="',dbulk,'"/> <!-- AKA envelope-->'));
fprintf(fid,'\t\t\t%s\n',strcat('<derode value="',derode,'"/> <!-- erosion distance; zero to disbale this step-->'));
fprintf(fid,'\t\t%s\n','</algorithm>');
fprintf(fid,'\t\t%s\n',strcat('<range file="',rrngFile,'"/>'));
fprintf(fid,'\t\t%s\n','<!--select ions for cluster core (aka solute)-->');
fprintf(fid,'\t\t%s\n','<core>');
fprintf(fid,'\t\t\t%s\n','<typelist>');
for i = 1:length(coreList)
    fprintf(fid,'\t\t\t%s\n',strcat('<atomtype symbol="',coreList{i},'"/>'));
end
fprintf(fid,'\t\t\t%s\n','</typelist>');
fprintf(fid,'\t\t%s\n','</core>');
fprintf(fid,'\t\t%s\n','<!--select ions for bulk (aka matrix)-->');
fprintf(fid,'\t\t%s\n','<bulk>');
fprintf(fid,'\t\t\t%s\n','<typelist>');
for i = 1:length(bulkList)
    fprintf(fid,'\t\t\t\t%s\n',strcat('<atomtype symbol="',bulkList{i},'"/>'));
end
fprintf(fid,'\t\t\t%s\n','</typelist>');
fprintf(fid,'\t\t%s\n','</bulk>');
fprintf(fid,'\t\t%s\n',strcat('<sizeclip nmin="',nmin,'"/>'));

%<!-- do we wish to keep unranged ions after analysing?-->
if unranged
    fprintf(fid,'\t\t%s\n','<unranged foroutput="true" forstats="true"/>');
end

if ischar(statFile)
    fprintf(fid,'\t\t%s\n',strcat('<clusterstats core="true" bulk="true" percluster="true" file="',statFile,'"/>'));
end
if(clustered)
    fprintf(fid,'\t\t%s\n',strcat('<clustered-pos file="',clustered,'" retain="true"/>'));
end
if(unclustered)
    fprintf(fid,'\t\t%s\n',strcat('<unclustered-pos file="',unclustered,'" retain="true"/>'));
end
if(clusterIndex)
    fprintf(fid,'\t\t%s\n',strcat('<clusterid file="',clusterIndex,'" offset="1"/>'));
end
fprintf(fid,'\t%s\n','</cluster>');
if ~strcmp(sweep,'none');
    fprintf(fid,'\t%s\n','</clustersweep>');
end
fprintf(fid,'%s\n','</posscript>');
fclose(fid);
disp(strcat('Wrote:',xmlName));
end
