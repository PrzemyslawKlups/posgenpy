% Runs multiple cluster searchs of random data and returns the cluster size
% distribution in number of core atoms in the core.
% A London Nov 2013
% number of runs to perform
% WARNING! Each run may take at least 20s to run!
tic
num_runs = 30; % number of random repeats
core = {'YO' 'TiO' 'Y' 'O'}; % core ions
fileId = 'someID_'; % something to identify the files
% file inputs
posFile = 'some.pos';
rrngFile = 'some.rrng';

% Set posgen dir
posDir = 'C:\Users\andy\APT\posgen\';
% copy rrng file to posgen dir
cpy(rrngFile,posDir);

[~,rrngName,rrngExt]= fileparts(rrngFile);
[~,posName,posExt]= fileparts(posFile);

dselection = 0.3:0.2:0.9;
xd = cell(size(dselection));
yd = xd;
for d = 1:length(dselection)
	% Run cluster search
	% Generate clusterOps.xml file:
    dmax = num2str(dselection(d));%'0.7';
    dbulk = '0'; % set to zero to skip, time saver!
    nmin = '3';
    statsFile = strcat(fileId,'clusteredStats.txt');
    disp('Generating posgen options xml file');
    posgenOptionGen('ops.xml',strcat(posName,posExt),strcat(rrngName,rrngExt),dmax,'1',dbulk,'0',nmin,core,statsFile,0,0,0,0);
    if(0) % for cluster search of real data
        % Run posgen
        disp('Running posgen');
        returncode = unix('posgen ops.xml');
        if returncode ~=0
            error('posgen died')
        end
        % get results (read statsFile)
        disp('Reading clustered stats file');
        [~, ~,~,~,~,~, clusteredN] = clusterStatsReader(statsFile, core);
    end
    xx = []; yy=[];
    unclusterResults = cell(1,num_runs);
    %figure; hold all
    for r = 1:num_runs
        % run search of random data
        statsFile = strcat(fileId,'unclusteredStats.txt');
        disp('Running posgen for unclustered search');
        posgenOptionGen('ops.xml',strcat(posName,posExt),strcat(rrngName,rrngExt),dmax,'1',dbulk,'0',nmin,core,statsFile,1,0,0,0);
        % Run posgen
        disp('Running posgen');
        returncode = unix('posgen ops.xml');
        if returncode ~=0
            error('posgen died')
        end
        % get results
        disp('Reading unclustered stats file');
        [~, ~,~,~,~,~, unclusterResults{r}] = clusterStatsReader(statsFile, core);
        ux = (str2double(nmin):(max(unclusterResults{r})+1));
        [uy] = histc(unclusterResults{r},ux-.5);
        %plot(ux,uy,'x');
        % store results
        xx = [xx ux];
        yy = [yy; uy];
    end
    xd{d} = xx;
    yd{d} = yy;
end
%set(gca,'yscale','log'); % **use log scale**
%xlabel('Core ions count');ylabel('Number of clusters');
%hold off
% plot with dmax
figure
loglog(xd{1},yd{1},'x')
hold all
for i = 2:length(dselection)
    loglog(xd{i},yd{i},'x')
end
hold off
% generate legend
[~,~,~,~,~,~,splitstring] = regexp(num2str(dselection,'%2.2f\t'),'\t');
legend(splitstring);
xlabel('Core ions count');ylabel('Number of clusters');
title('Cluster Size Distribution for various D_{MAX} nm')
toc