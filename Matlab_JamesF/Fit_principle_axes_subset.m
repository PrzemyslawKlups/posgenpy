function [allaxes, sizenatoms, clustercenter,P] = Fit_principle_axes_single_cluster(indxClrPos,fignum,COIID,plotelip) 
%Takes an indexed cluster posfile, runs principal component analysis, to
%determine principal direction of the cluster axes, and their lengths.
%plots a 3D visualisation of the fit and outputs the numerical data%
%Inputs%
%indexClrPos - indexed cluster pos file%
%COIID - Cluster of interest ID%
%OutPuts%
%allaxes - allaxes(nthcluster, xyzl, axes 123) - matrix of princpal axes for each cluster and it's length (diameter not radius)%
%axis 1 > axis 2 > axis 3%

[x, y, z, m, nb] = readpos(indxClrPos);
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
figure (fignum)
clf
nthcluster = 0;
allaxes = zeros(size(COIID,2), 4, 3);
sizenatoms = zeros(size(COIID,2),1);
clustercenter = zeros(size(COIID,2),3);
for clusterID = COIID
    nthcluster = nthcluster + 1;
X = 0;
Y = 0;
Z = 0;
natomsincluster =0;
for i = 1:size(pos,1)
    if pos(i,4) == clusterID
        X = [X; pos(i,1)];  
        Y = [Y; pos(i,2)];
        Z = [Z; pos(i,3)];
        natomsincluster = natomsincluster + 1;
    end
end
sizenatoms(nthcluster) = natomsincluster; 
X = X(2:size(X));
Y = Y(2:size(Y));
Z = Z(2:size(Z));
%following principal component analysis from wikipeadia Feb 18, some
%crossing of nomenclature%
P = [X, Y, Z]; %1 one in p for each ion in the cluster%
c = mean(P)'; %column vector of centre of mass of cluster (not weighted by ion mass)%
clustercenter(nthcluster,:) = c; %output matrix of all cluster centre points
%egien vlaue method of finding principla components, no need to do this however as Matlab has a principle componet function%    
    %h = ones(natomsincluster,1);
    %B = P - h*c'; %subtract centre from each point, so COM is now the origin%
    %C = (B'*B)/(natomsincluster -1); %Covariance matrix%%Minus 1 comes from Bessels correction, don't know what it is%
    %[evec,eval]= eig(C); %eigne vector of the coraince matrix are the principla axes of the cluster. 
[axes,score,latent,tsquared]= pca(P); %axes are the priciple axes of the cluster,
%score is something I don't understand, latent is the variance along each
%of the principle axes,tsquared may relate to getting a measure of goodness of fit%
lengths = 2*latent.^0.5;
%4 x times standard diavation should encompass ~95% of the points (so long
%as they are normally distributed?)
%NB at this stage lengths refers to radii
%tsquared may relate to getting a measure of goodness of fit%
%[evt,eval] = eig((ccpos(:,1:3)'*ccpos(:,1:3))/size(ccpos,1));
%later plotting needs, c centre of ellipse, A, matrix satifying
%(v-c)'*A*(v-c) = 1 and lengths of the vectors.
A = axes*diag(lengths.^-2)*inv(axes);
lengths = diag(lengths);

%draw fit
if strcmp(plotelip,'Yes') == 1
    scatter3(X,Y,Z,10,'.b')%plots the points in 3D, like IVAS or 3Dpict%
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    if natomsincluster < 1000
        nsteps = 100;
    else
        nsteps = 50;
    end
    [x,y, z] = meshgrid( linspace(c(1)-max(abs(X)),c(1)+max(abs(X)),nsteps), linspace(c(2)-max(abs(Y)),c(2)+max(abs(Y)),nsteps),linspace(c(3)-max(abs(Z)),c(3)+max(abs(Z)),nsteps));
    Ellipsoid = zeros(nsteps, nsteps, nsteps);
    for i = 1:nsteps
        for j = 1:nsteps
            for k = 1:nsteps
                v = [x(i,j,k); y(i,j,k); z(i, j, k)]; %v is column vector of [x; y; z] %
                Ellipsoid (i,j,k) = (v - c)' * A * (v-c); 
            end
        end
    end
    hold on
    p = patch( isosurface( x, y, z, Ellipsoid, 1) ); %The Ellipse statifises (v-c)'*A*(v-c) = 1 so plot on isosurface at the value of 1%
    hold off;
    set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
    axis vis3d equal;
    camlight;
    alpha 0.3 %ellipse transparency
    lighting phong;
end
hold on
for i = 1:3
    principaxis = [c - axes(:,i).*lengths(i,i), c + axes(:,i).*lengths(i,i)];
plot3(principaxis(1,:),principaxis(2,:),principaxis(3,:))
end
hold off
allaxes(nthcluster,1:3,:) = axes; %all the axes directions%
for i = 1:3
allaxes(nthcluster,4,i) = 2*lengths (i,i);%NB output is the diameter now while lengths is actually the radius%
end
end
axis vis3d equal;