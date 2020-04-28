function [z_Plot,x_dmax,y_nmin,rawData] = posgenDmaxNminSweep(rangeFile, posFile, dmax, coreIons, relabel, proot)
%[z_Plot,x_dmax,y_nmin,rawData] = posgenDmaxNminSweep(rangeFile, posFile, dmax, coreIons, proot)
% rangeFile = rrngFile
% dmax, vector of dmaxs to use, requires constant difference, .1  .2 .3 etc
% coreIons = cell array of natural names {'Y','YO','Cu'} etc.
% relabel = true/false should the masses be relabelled?
% proot = the root directory to the posgen exe
% A London Sept 2016

if nargin < 6
    proot = 'C:\cygwin64\home\andy\posgen\';
end
timeStr = datestr(now,30);
xmlName = [proot timeStr '_sweepTest.xml'];
statFile = [timeStr '_sweepTest_stats.txt'];

[coreList, bulkList]=posgenOptionGen2(...
    xmlName,posFile,rangeFile,dmax,'1','0','0','2',coreIons,...
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

rawData = cell(length(filesToSearch),2);

for i = 1:length(filesToSearch)
    sfile = [filesToSearch{i} '/' statFile];
    % find dmax value
    dmaxValue = str2double(regexprep(filesToSearch{i},'sweep-dmax-',''));
    [elements1, countsDecomp1,counts1,position1,ionType1,rg1, nTable1] = clusterStatsReader( sfile, coreIons );
    rawData{i,1} = dmaxValue;
    rawData{i,2} = nTable1;
end
allN = vertcat(rawData{:,2});
[~,edges] = histcounts(allN,'binmethod','integers');
rawData2 = zeros(length(filesToSearch),length(edges)-1);
for i = 1:length(filesToSearch)
    [rawData2(i,:)] = histcounts(rawData{i,2},edges);
end
rawData3 = fliplr(cumsum(fliplr(rawData2),2)); % how many clusters with at least nmin ions

z_Plot = rawData3;
x_dmax = reshape(cell2mat(rawData(:,1)),1,[]);
y_nmin = reshape(min(allN):max(allN),1,[]);
% make a nice plot
%[X,Y] = meshgrid(reshape(x_dmax,1,[]), reshape(y_nmin,1,[]));
%surf(Y,X,min(a,300)','FaceColor','interp','EdgeColor','none','FaceLighting','gouraud')