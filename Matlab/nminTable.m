function [nTable] = nminTable(statsFile,coreIons)
% [nTable] = nminTable(statsFile,coreIons)
% Reads a cluster-stats file from posgen and returns a list of the total
% number of core atoms in each cluster as a 1D list
% TODO: Merge with clusterStatsReader
% A. London May 2013

% Read statsFile

[elements,countsDecomp,counts,position,ions,rg] = clusterStatsReader(statsFile,coreIons);

% alloc a results table
if(~isempty(data))
    nTable = zeros(length(data),1);
else
    error(strcat('Failed to import data file:',statsFile));
end
for c=1:length(ions) % loop over columns checking for coreAtoms
    %starts at 4, cause the first 3 cols are x,y,z
    for i=1:length(coreIons) % loop over each core ion
        if(strcmp(coreIons{i},ions{c}))
            % found a match, add to total
            nTable = nTable + counts(:,c);
        end
    end
end

end