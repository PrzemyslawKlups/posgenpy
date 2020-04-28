function [compSummary, ionType] = posgenErosionCompSweep(rangeFile, posFile, dmax, nmin, dbulk, derosion, coreIons, relabel, proot)
%[compSummary, ionType] = posgenErosionCompSweep(rangeFile, posFile, dmax, nmin, dbulk, derosion, coreIons, relabel, proot)
% rangeFile = rrngFile
% dmax = fixed value of maximum separation
% nmin = minimum number of core ions to be clustered
% dbulk = envoloping distance
% derosion = vector of erosion distances (must be constant step)
% coreIons = cell array of natural names {'Y','YO','Cu'} etc.
% relabel = true/false should the masses be relabelled?
% proot = the root directory to the posgen exe
% returns a matrix which is derode long and ion types wide
% A London Sept 2016

if nargin < 9
    proot = 'C:\cygwin64\home\andy\posgen\';
end
timeStr = datestr(now,30);
xmlName = [proot timeStr '_sweepTest.xml'];
statFile = [timeStr '_sweepTest_stats.txt'];

[coreList, bulkList]=posgenOptionGen2(...
    xmlName,posFile,rangeFile,dmax,'1', dbulk, derosion, nmin, coreIons,...
    statFile,relabel,0,0,0,false);
% delete existing sweep folders
files = dir('sweep*');
dirs = cell2mat({files.isdir})==1;
names = {files.name}';
filesToDelete = names(dirs);
for i = 1:length(filesToDelete)
    movefile(filesToDelete{i},['oldSweeps/' filesToDelete{i}]);
end
disp('Running posgen:');
returncode = unix([proot 'posgen.exe ' xmlName]);
if returncode ~=0
    error('posgen died')
end
% extract statistics from folders
% find folders:
files = dir('sweep*');
dirs = cell2mat({files.isdir})==1;
names = {files.name}';
filesToSearch = names(dirs);

if isempty(filesToSearch)
    % only a single file was created
    % read this instead
    % has the name: statFile
    [elements, countsDecomp1,counts1,position1,ionType,rg1, nTable1] = clusterStatsReader( statFile, coreList );
    rawData{1,1} = '';
    rawData{1,2} = counts1;
else
    
    rawData = cell(length(filesToSearch),2);
    
    for i = 1:length(filesToSearch)
        sfile = [filesToSearch{i} '/' statFile];
        % find dmax value
        dmaxValue = str2double(regexprep(filesToSearch{i},'sweep-derode-',''));
        [elements, countsDecomp1,counts1,position1,ionType,rg1, nTable1] = clusterStatsReader( sfile, coreList );
        rawData{i,1} = dmaxValue;
        rawData{i,2} = counts1;
    end
end
allComp = cellfun(@(x) sum(x,1), rawData(:,2),'uni',0);
countSummary = vertcat(allComp{:});
compSummary = bsxfun(@times, countSummary, 1./sum(countSummary,2));