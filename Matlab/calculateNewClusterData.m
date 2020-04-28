function [newcsv] = calculateNewClusterData(newFile,rrngFile,clrPos,indxClrPos,matPos,clrSettings, clusterIDoffset)
%[newcsv] = calculateNewClusterData(newFileName,rrngFile,clrPos,indxClrPos,matrixPosFile,clrSettings, clusterIDoffset)
% Makes a new cluster data file (.csv) similar to that produced by IVAS
% using a rrng file, cluster.pos file, cluster.indexed.pos file, the pos 
% file of the matrix (everything not clustered) and a
% definition of the cluster settings, a cell formatted like so:
%{ROI;core;dmax;order;nmin;dbulk;derode}
% You can generate this using: clusterSettingsNumeric
%[dmax,dbulk,derode,nmin,core,count,numberofions,order,ROI] = ... 
% clusterSettingsNumeric('C:\Users\andy\APT\clusterAnalysis\R14_16469-v01_YYOTiO_80.csv');
%clrSettings = {ROI;core;dmax;order;nmin;dbulk;derode};
% clusterIDoffset = the offset to add to the loaded cluster IDs, 0 for IVAS
%   1 for posgen.
% This returns a cell formatted like the csv (and just as hard to read!) as
% well as writing a new csv file called 'newFile'
% A London Nov 2014

% clusterID offset
%clusterIDoffset = 1; % 1 for posgen, 0 for ivas
[element_num, range_num, elements, ranges] = rangeReader(rrngFile);
[ix,iy,iz,clusterIds] = readpos(indxClrPos);
iPos = [ix' iy' iz']; 
clear ix iy iz
clusterIds = clusterIds + clusterIDoffset;
totalClusters = round(max(clusterIds)); % how many clusters?

%% Range Data in all, ranged and solute/core
% all
[x,y,z,m] = readpos(clrPos);
aPos = [x' y' z' m']; % all data inc un-ranged
% register aPos and indexed (cluster IDs)
[~,indices]=ismemberNoNAN(aPos(:,1:3),iPos,'rows');
clusterIds = clusterIds(indices);
clear x y z iPos
% core/solute - only masses that match core ranges
% which ranges are core ranges?
coreRange = ranges(ismember(cell2mat(ranges(:,5:end)),ions2ionTable(ionStr2ions(clrSettings{2}),elements),'rows'),:);
rangeTableCore = cell2mat(coreRange(:,[1 2 5:end]));
[ionNumCore,~] = rangeM(m,rangeTableCore,elements);
cPos = aPos(ionNumCore>0,:); % core-ion only clusters
cIDs= clusterIds(ionNumCore>0);
clear ionNumCore
% ranged
rangeTable = cell2mat(ranges(:,[1 2 5:end])); %type double 
[ionNum,~] = rangeM(m,rangeTable,elements);
rPos = aPos(ionNum>0,:); % ranged-ion only clusters
rIDs= clusterIds(ionNum>0);

%% Calculate new cluster stats
% solute
% ionic range table:
[ionicRanges,ions] = ele2ionicRanges(rangeTable,elements); % used later for composition
numOfIons = length(ions);
clusterData = cell(totalClusters,71+(numOfIons+1)*3);

for c = 1:totalClusters
    for i=1:3 % one each for 1Solute/2Ranged/3Total
        co = [5 17 29]; % column offset value
        co = co(i); % this is so you only have to change the offsets here ^
        if i == 1 % solute/core basis
            ccpos = cPos(cIDs==c,:); % current cluster pos [x y z m]
            clusterData{c,2} = size(ccpos,1); % set total solute ions
        elseif i == 2 % ranged basis
            ccpos = rPos(rIDs==c,:); % current cluster pos [x y z m]
            clusterData{c,3} = size(ccpos,1); % set total ranged ions
        elseif i == 3 % all basis
            ccpos = aPos(clusterIds==c,:); % current cluster pos [x y z m]
            clusterData{c,4} = size(ccpos,1); % set total ions
        end
        clusterData{c,1+co} = sum(ccpos(:,1).*ccpos(:,4))./sum(ccpos(:,4));
        clusterData{c,2+co} = sum(ccpos(:,2).*ccpos(:,4))./sum(ccpos(:,4));
        clusterData{c,3+co} = sum(ccpos(:,3).*ccpos(:,4))./sum(ccpos(:,4));
        [~,rx,ry,rz]=radiusGyration(ccpos(:,1),ccpos(:,2),ccpos(:,3),ccpos(:,4));
        clusterData{c,4+co} = rx;
        clusterData{c,5+co} = ry;
        clusterData{c,6+co} = rz; % radii solute nm
        clusterData{c,7+co} = (4/3)*pi*rx*ry*rz; % V_Rg nm^3
        cmin = min(ccpos(:,1:3),[],1);
        cmax = max(ccpos(:,1:3),[],1);
        clusterData(c,(8:10)+co) = num2cell(mean([cell2mat(clusterData(c,(1:3)+co))-cmin; cmax-cell2mat(clusterData(c,(1:3)+co))],1)); % Extent X,Y, Z
        clusterData{c,11+co} = prod(cell2mat(clusterData(c,(8:10)+co)))*pi*(4/3); % V_Extent nm^3
        
        %co = co - 3; % column adjustment
        % transform to the principal coords of the cluster
        [evt,eval] = eig((ccpos(:,1:3)'*ccpos(:,1:3))/size(ccpos,1));
        ctpos = [ccpos(:,1:3)*evt ccpos(:,4)]; % current transformed pos
        clusterData{c,37+co} = clusterData{c,1+co};
        clusterData{c,38+co} = clusterData{c,2+co};
        clusterData{c,39+co} = clusterData{c,3+co};
        [~,rx,ry,rz]=radiusGyration(ctpos(:,1),ctpos(:,2),ctpos(:,3),ctpos(:,4));
        clusterData{c,40+co} = rx;
        clusterData{c,41+co} = ry;
        clusterData{c,42+co} = rz; % radii solute nm
        clusterData{c,43+co} = (4/3)*pi*rx*ry*rz; % V_Rg nm^3
        clusterData{c,44+co} = (max(ctpos(:,1))-min(ctpos(:,1)))/2; % Extent X
        clusterData{c,45+co} = (max(ctpos(:,2))-min(ctpos(:,2)))/2; % Extent Y
        clusterData{c,46+co} = (max(ctpos(:,3))-min(ctpos(:,3)))/2; % Extent Z
        clusterData{c,47+co} = prod(cell2mat(clusterData(c,(44:46)+co)))*pi*(4/3); % V_Extent nm^3
    end
    % compositions
    [~,ionicCounts] = massQuant(ccpos(:,4),ionicRanges);
    clusterData(c,78:(length(ions)+77)) = cellfun(@(x) sprintf('%6.4f%%',x), num2cell(100*ionicCounts./clusterData{c,3}),'uni',0); % range%
    clusterData(c,(length(ions)+77+2):(2*length(ions)+77+1)) = cellfun(@(x) sprintf('%6.4f%%',x), num2cell(100*ionicCounts./clusterData{c,4}),'uni',0); % total%
    clusterData(c,(2*length(ions)+77+3):(3*length(ions)+77+2)) = num2cell(ionicCounts); %ionic counts
end
%% Make matrix data
% All matrix data
[matrixMasses] = loadMasses(matPos);
totalMatIons = length(matrixMasses);
% Core/solute ions only (count)
[ionNumMat,~] = rangeM(matrixMasses,rangeTableCore,elements);
totalMatSoluteIons = sum(ionNumMat>0);
clear ionNumMat
% Ranged data, mass list and count
[ionNum,~] = rangeM(matrixMasses,rangeTable,elements);
rangedMatrix = matrixMasses(ionNum>0);
totalMatRangedIons = length(rangedMatrix);
clear ionNum
% ionic counts
[~,ionicCounts] = massQuant(matrixMasses,ionicRanges);
% matrix table (cell)
mn = cell(1,size(clusterData,2));
mn{2} = totalMatSoluteIons;
mn{3} = totalMatRangedIons;
mn{4} = totalMatIons;
% 78 is the column number where the comp data starts
mn(78:(length(ions)+77)) =                      cellfun(@(x) sprintf('%6.4f%%',x), num2cell(100.*ionicCounts./totalMatRangedIons) ,'uni',0);
mn((length(ions)+77+2):(2*length(ions)+77+1)) = cellfun(@(x) sprintf('%6.4f%%',x), num2cell(100.*ionicCounts./totalMatIons) ,'uni',0);
mn((2*length(ions)+77+3):(3*length(ions)+77+2))=num2cell(ionicCounts);

%% make new csv style file
% header
newcsv(1,1:2) = ['ROI' clrSettings(1)];
newcsv(2,1:(1+size(clrSettings{2},2))) = ['Ion(s)' clrSettings{2}];
newcsv(3,1:2) = ['d-max (nm)' clrSettings(3)];
newcsv(4,1:2) = ['Order (ions)' clrSettings(4)];
newcsv(5,1:2) = ['N-min (ions)' clrSettings(5)];
newcsv(6,1:2) = ['L (nm)' clrSettings(6)];
newcsv(7,1:2) = ['d-erosion (nm)' clrSettings(7)];
newcsv(9,1:2) = {'Cluster Count',totalClusters};
%% Header column titles
% total counts
h = cell(1,size(clusterData,2));
h(2:4) = strcat({'Solute','Ranged','Total'},{' Ions'});
% shape descriptors
shapeDes = {'Center_x (nm)','Center_y (nm)','Center_z (nm)','R_gx (nm)','R_gy (nm)','R_gz (nm)','V_Rg (nm^3)','Extent_x (nm)','Extent_y (nm)','Extent_z (nm)','V_Extent (nm^3)'};
cstart = 6;
for i = 1:3
    if i==1
        a = ' Solute';
    elseif i==2
        a = ' Ranged';
    elseif i==3
        a = ' Total';
    end
    for j = 1:2
        if j==1
            b = ''; % add nothing for original data
        elseif j==2
            b = ''''; % add ' for transformed data
        end
        
        hstart = cstart + (length(shapeDes)+1)*(i-1) + (length(shapeDes)+1)*3*(j-1);
        hend   = cstart + (length(shapeDes)+1)*(i) + (length(shapeDes)+1)*3*(j-1) - 2;
        
        h(hstart:hend) = strcat(shapeDes,a,b);
        
    end
end
% composition data headers
h(78:(length(ions)+77)) = strcat(ions',' % Ranged');
h((length(ions)+77+2):(2*length(ions)+77+1)) = strcat(ions',' % Total');
h((2*length(ions)+77+3):(3*length(ions)+77+2))= strcat(ions',' Count');
%% headers
newcsv(11,1:size(h,2)) = h;
%% matrix
newcsv(12,1:(size(h,2))) = ['Matrix' mn(2:end)];
%% cluster names
newcsv(13:(totalClusters+12),1) = cellfun(@(x) sprintf('Cluster %d',x),num2cell(1:totalClusters)','uni',0);
% % cluster data
newcsv(13:(totalClusters+12),2:(size(h,2))) = clusterData(:,2:end);
% % write new cluster data file
disp(['Writing: ' newFile]);
cell2csv(newFile,newcsv);
end