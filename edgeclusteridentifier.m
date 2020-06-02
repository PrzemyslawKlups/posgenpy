function [centralCOIID,edgeCOIIDsmall,edgeCOIIDlarge,V] = edgeclusteridentifier(indxClrPos,wholeDataPos,Dthreshold,sizethreshold,fignum)  

%This function indentifies clusters which are on the edge of a datasets by
%creating an alpha hull, more details about that below. It plots 6 graphs,
%showing the hull, atom maps of which clusters are categoried as what, axes
%fit for each category and then a morphology graph.

%Does take time%

%Inputs%
%indxClrPos- indexed cluster posfile

%wholeDataPos - whole reconstruction from which the indexed clusters come
%from

%Dthreshold -maximum distance that any point in the cluster can be away
%from one of the vertices in the alphahull, in nm. ,default 0.5 nm

%sizethreshold - idea is that small clusters on the edge aren't as
%interesting, high amount of error, but large one, i.e. needles that hit
%the edge and have definite shape and orientation are, but I still need to
%clarify they're on the edge. default is 1000 detected atoms

%V - volume of whole dataset, can use for calculating number densities

%Outpus%

%Comments%
%Colour is hard coded, could be changed but this is likely an internal
%functin so not worth it.
%Alpha hull part is slow, it's a computationally expensive job to display
%as well as make the hull. There is a hard coded sampling factor.
%it would be quicker to use a pre-determined value of alpha but I've gone
%with optimising alpha as this seems more generalised

if strcmp(Dthreshold,'default') ==1
 Dthreshold = 0.5;
end
if strcmp(sizethreshold,'default') ==1
 sizethreshold = 1000;
end
[xw, yw, zw, mw, nbw] = readpos(wholeDataPos);
wholepos = [xw; yw; zw]'; %wholepos only has three columns, I've left out mass to charge as this is currently unimportant%
wholepos = double(wholepos);
%samplewpos = datasample(wholepos,ceil(size(wholepos,1))); %sample the dataset to reduce size by factor of 10,000 to speed up computation%
samplingfactor = 10;
nsamples = ceil(size(wholepos,1)/samplingfactor) -1; %number of samples 
samplerows = [1:1:nsamples]*samplingfactor; %which rows are being sampled, sampling is evenly spaced not random%
samplewpos = zeros(nsamples,3);
for i = 1:nsamples
    
    samplewpos(i,:)= wholepos(samplerows(i),:);
end
tipsurface = alphaShape(samplewpos,'HoleThreshold',10^100);
figure(fignum)
clf
plot(tipsurface)
V = volume(tipsurface);

[x, y, z, m, nb] = readpos(indxClrPos);
clusterpos = [x; y; z; m]';
clusterpos = double(clusterpos);
numberclusters = 0;
clusterIDS = 0;
for i = 1:max(clusterpos(:,4))
    if ismember(i,clusterpos(:,4)) == 1
       numberclusters = numberclusters + 1;
       clusterIDS = [clusterIDS , i];
    end
end
clusterIDS = clusterIDS (2:numberclusters+1);
figure (fignum + 1)
clf
nthcluster = 0;
allaxes = zeros(numberclusters, 4, 3);
sizenatoms = zeros(numberclusters,1);
clustercenter = zeros(numberclusters,3);
edgeCOIIDsmall = 0;
edgeCOIIDlarge = 0;
centralCOIID = 0;
for clusterID = clusterIDS
    nthcluster = nthcluster + 1;
X = 0;
Y = 0;
Z = 0;
natomsincluster =0;
for i = 1:size(clusterpos,1)
    if clusterpos(i,4) == clusterID
        X = [X; clusterpos(i,1)];  
        Y = [Y; clusterpos(i,2)];
        Z = [Z; clusterpos(i,3)];
        natomsincluster = natomsincluster + 1;
    end
end
sizenatoms(nthcluster) = natomsincluster; 
X = X(2:size(X));
Y = Y(2:size(Y));
Z = Z(2:size(Z));
%scatter3(X,Y,Z)%would plot the points in 3D, like IVAS or 3Dpict%
xlabel('X')
ylabel('Y')
zlabel('Z')
%following principal component analysis from wikipeadia Feb 18, some
%crossing of nomenclature%
P = [X, Y, Z]; %1 one in p for each ion in the cluster%

[I,D] = nearestNeighbor(tipsurface,P); %indices of nearest point in tipsurface for each point in P, and distance D%
hold on
if any(D < Dthreshold)
    %cluster is close to the edge%
    if natomsincluster < sizethreshold
        edgeCOIIDsmall = [edgeCOIIDsmall, clusterID];
        scatter3(X,Y,Z,'r.')
    else
        edgeCOIIDlarge = [edgeCOIIDlarge, clusterID];
        scatter3(X,Y,Z,'y.')
    end
else
    centralCOIID = [centralCOIID, clusterID];
    scatter3(X,Y,Z,'b.')
end
end
edgeCOIIDsmall = edgeCOIIDsmall(2:size(edgeCOIIDsmall,2));
edgeCOIIDlarge = edgeCOIIDlarge(2:size(edgeCOIIDlarge,2));
centralCOIID = centralCOIID(2:size(centralCOIID,2));
axis vis3d equal;

% [allaxessmalledge, sizenatomssmalledge, clustercenter,P] = Fit_principle_axes_single_cluster(indxClrPos,fignum+2,edgeCOIIDsmall,'No');
% [allaxeslargeedge, sizenatomslargeedge, clustercenter,P] = Fit_principle_axes_single_cluster(indxClrPos,fignum+3,edgeCOIIDlarge,'No');
% [allaxescentral, sizenatomscentral, clustercenter,P] = Fit_principle_axes_single_cluster(indxClrPos,fignum+4,centralCOIID,'No');
% plotdatalargeedge =  morphologygraph(allaxeslargeedge,sizenatomslargeedge,fignum+5,'y',1,'No','Yes');
% plotdatasmalledge =  morphologygraph(allaxessmalledge,sizenatomssmalledge,fignum+5,'r',1,'No','No');
% plotdatacentral =  morphologygraph(allaxescentral,sizenatomscentral,fignum+5,'b',1,'atoms','No');
% legend('edgelarge','edgesmall','central','Location','eastoutside')