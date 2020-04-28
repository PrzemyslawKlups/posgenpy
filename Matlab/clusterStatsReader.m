function [elements, countsDecomp,counts,position,ionType,rg, nTable] = clusterStatsReader( statsFile, coreIons )
%[elements, countsDecomp,counts,position,ionType,rg, nTable] = clusterStatsReader( clusterFile, coreIons )
%CLUSTERSTATSREADER This reads the clusterStats.txt file produced by posgen
% returns the array of elements (column headings, cell array) and
% the decomposed counts (matrix)

% import data
data=dlmread(statsFile,'\t',1,0); % open file skipping header
fid=fopen(statsFile);
tline = fgets(fid);
fclose(fid);
textdata = strsplit(tline,'\t','DelimiterType','RegularExpression');
% textdata = headers
% data = numbers


ionType = textdata(1,4:(end-1));% does not ignore un-ranged
counts = data(:,4:(end-1));     %counts = cell2mat(counts);
position = data(:,1:3);         %position = cell2mat(position);
rg = data(:,end);               %rg = cell2mat(rg);

% Report num_clusters:
disp(strcat('Found:',num2str(length(data(:,1))),' clusters'));

% For decomposition
ions = ionStr2ions(ionType(1:end-1)); % don't pass in unranged
[ionTable,elements] = ions2ionTable(ions);

% decompose counts by doing matrix muliplication
countsDecomp = counts(:,1:end-1)*ionTable;

nTable = zeros(length(counts(:,1)),1);
% Make a column of core ions per cluster:
for c=1:length(ionType) % loop over columns checking for coreAtoms
    %starts at 4, cause the first 3 cols are x,y,z
    for i=1:length(coreIons) % loop over each core ion
        if(strcmp(coreIons{i},ionType{c}))
            % found a match, add to total
            nTable = nTable + counts(:,c);
        end
    end
end

% column headings is just 'elements'
end