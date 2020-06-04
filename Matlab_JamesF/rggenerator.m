function [rg,COIID] = rggenerator(indxClrPos,COIID)
%this will output a rg based on all ions in the cluster

%Inputs - indexed cluster posfile
%COIID - cluster of interest ID, can be set to 'all'

%copied the fit pricimple axes code and editted%
[x, y, z, m, ~] = readpos(indxClrPos);
if strcmp (COIID,'all') == 1
    COIID = unique(m);
end  
pos = [x; y; z; m]';
pos = double(pos);
numberclusters = 0;
n = 0;
for i = 1:max(pos(:,4))
    if ismember(i,pos(:,4)) == 1
       numberclusters = numberclusters + 1;
       n = [n , i];
    end
end
nthcluster = 0;
rg = zeros(size(COIID,2),1);
for clusterID = COIID
    nthcluster = nthcluster + 1;
X = 0;
Y = 0;
Z = 0;
for i = 1:size(pos,1)
    if pos(i,4) == clusterID
        X = [X; pos(i,1)];  
        Y = [Y; pos(i,2)];
        Z = [Z; pos(i,3)];
    end
end 
X = X(2:size(X));
Y = Y(2:size(Y));
Z = Z(2:size(Z));
P = [X, Y, Z]; %1 one in p for each ion in the cluster%
c = mean(P); %row vector of centre of mass of cluster (not weighted by ion mass)%
r = rssq(P-c,2);%distance from centre of mass%
rg(nthcluster,:) = rms(r); %output matrix of all cluster's radii of gyration
end