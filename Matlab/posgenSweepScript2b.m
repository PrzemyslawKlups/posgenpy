% posgen sweep script to produce composition 3D plot
% 	load saved data for example: load('posgenSweepScript2b_results.mat')

% firstly automate the dmax-nmin selection (using a dmax sweep)
posFile = 'Dataset2.pos';
rangeFile = 'Dataset1.rrng';
%dmax = [0.3:0.01:0.7];
coreIons = {'Cu'};

% do a sweep clustering search then extract the average composition as a
% function of erosion and includes distance. The includes distance is set
% below and the erosion will be from 0 to the current value of includes
includes = 0.01:0.01:0.25;
X = cell(length(includes),1); % includes value
derosion = cell(length(includes),1);
result = cell(length(includes),1);
for i = 1:length(includes)
    tic
    disp(['*****' num2str(i) '/ ' num2str(length(includes))])
    derosion{i} = 0:0.01:includes(i);
    [compSummary,ionTypes] = posgenErosionCompSweep(rangeFile, posFile, '1.0', '50', num2str(includes(i)), derosion{i}, coreIons, false);
    result{i} = compSummary(:,2);
    X{i} = ones(size(compSummary(:,2)))*includes(i);
    toc
end
X = vertcat(X{:}); % includes value
Y = horzcat(derosion{:})'; % erosion value
Z = vertcat(result{:}); % average composition of clusters

%% 3D scatter plot: includes-erosion-composition
scatter3(X,Y,Z);
xlabel('Includes dist (nm)');
ylabel('Erosion dist (nm)');
zlabel('Mean cluster comp. (at%)');
title('Dmax = 1.0, Nmin = 50');

%% 3D surface plot: includes-erosion-composition
figure
tri = delaunay(X,Y);
h = trisurf(tri, X, Y, Z);
xlabel('Includes dist (nm)');
ylabel('Erosion dist (nm)');
zlabel('Mean cluster comp. (at%)');
title('Dmax = 1.0, Nmin = 50');

%% coloured surface plot of 2nd diff gradient
% (not sure if this is useful)
[ux,~,ix] = uniquetol(X);
[uy,~,iy] = uniquetol(Y);
idx = sub2ind([max(ix) max(iy)],ix,iy);
F = nan([max(ix) max(iy)]);
F(idx) = Z;

figure
surf(F,gradient(gradient(F))*50)
caxis([-0.1 0.1])