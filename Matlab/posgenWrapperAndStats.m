% MATLAB Wrapper for using posgen.exe using posgenOptionGen.m
posgenDir = 'C:\cygwin64\home\andy\posgen\';
posRoot = 'F:\APT\R83_00675\recons\recon-v01\default\';
xmlName = [posgenDir 'inOps.xml']; % full path.name.ext
posFile = [posRoot 'R83_00675-v01.pos']; % full path.name.ext
rrngFile = [posRoot 'R83_00675-v01_new.cluster.rrng']; % full path.name.ext
dmax = '1.1';
knn = '1';
dbulk = dmax;
derode= '0.5';
nmin = '10';
core = {'YO' 'Y'};
statFile = [posRoot 'R83_00675-v01_new2_stats.csv'];
relabel = 0;
clustered = [posRoot 'R83_00675-v01_new2.cluster.pos'];
unclustered = [posRoot 'R83_00675-v01_new2.matrix.pos'];
clusterIndex = [posRoot 'R83_00675-v01_new2.cluster.indexed.pos'];

posgenOptionGen(xmlName,posFile,rrngFile,dmax,knn,dbulk,derode,nmin,core,statFile,relabel,clustered,unclustered,clusterIndex)

oldDir = cd(posgenDir); % cd to posgen directory (full path)

disp('Running posgen:');
returncode = unix([posgenDir 'posgen.exe ' xmlName]);
if returncode ~=0
    error('posgen died')
end

cd(oldDir); % change directories back

% make IVAS like cluster stats csv file
clrSettings = {'Top Level ROI-posgen';core;dmax;knn;nmin;dbulk;derode};
[newcsv] = calculateNewClusterData([posRoot 'R83_00675-v01_new2.csv'],...
    rrngFile,[clustered],[clusterIndex],[unclustered],clrSettings);