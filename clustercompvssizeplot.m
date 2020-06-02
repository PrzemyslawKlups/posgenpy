function [plotdata,averagecomp,clustersize] = clustercompvssizeplot(clusterstatsfile,numerIon,denomIon,scalespots,clearfig,fignum,nbins,linearscale,xvalue,Colour,erbar,plotmean)
%clustercompvssizeplot -plots composition or ratio of elements in
%clusters/precipiates as a function of size

%Inputs%
%clusterstatsfile - clusterstats file ouputted by posgen- NEEDS TO BE
%ALREADY FILTERED BY NMIN ALREADY
%numerIon - numerator ion.
%denomIon - denominator ion/ions for which you use to calculate composition
%or ratio - set as 'one' to allow plotting of just a number of ions or
%r_gyration. If doing this you and want to use nions as xvalue, may need to generate this using this function beforehand 
%scalespots - option to scale spot size by something, such as number of
%clusters in that bin
%clearfig - clear figure 'Yes' or 'No' - use 'No' to plot multiple datasets
%on the same graph.
%fignum - number of figure plotting on
%nbins - number of bins, or vector of binends to use. If not supplied the
%default is log10 spaced bins.
%Colour - spot colour
%linearscale - linear factor by which to scale spot size
%xvalue - rg for radius of gyration, nions, for numer of denumrator ions,
%ratio for sum of numerator and denomiator ions, 'Z' for tip axis Z value or a column vector of
%different values%
%erbar - currently Yes gives error bars of std error. Could maybe make it
%an optional std dev or std er.

%Known problems
%If doing a ratio and a cluster has none of the denominator ion in it, this
%will make that entire bin's average ratio Inf. 

%Example inputs
% clusterstatsfile = 'C:\Users\James\Documents\MATLAB\posgen\2ndyear_talk\v33_iPA_7025\clusterstats-dmax-0.7.txt';
% numerIon = {'Mg'};
% denomIon = {'Mg', 'Si', 'Cu'};
% scalespots = 'Yes';
% clearfig = 'Yes';
% fignum = 1;
% nbins = 20;
% Colour = 'r';
% linearscale = 1;
% xvalue = 'nions';

% import data - section taken from Andy's read clusterstats function
fid = fopen(clusterstatsfile,'r'); % open file
% Read it a line at a time
tline = fgets(fid);
% Replace tabs with commas
newStr = regexprep(tline, '\t', ',');
newStr = regexprep(newStr, '\n', ',');
% Check for duplicates
newStr = regexprep(newStr, ',,', ',');
C = textscan(newStr,'%s','Delimiter',',');
textdata = C{1}';
data = [];
% read the numeric data
while ischar(tline)
    % replace tabs with commas
    newStr = regexprep(tline, '\t', ',');
    newStr = regexprep(newStr, '\r', '');
    newStr = regexprep(newStr, '\n', '');
    % check for duplicates
    newStr = regexprep(newStr, ',,', ',');
    % Translate into numeric
    C = textscan(newStr,'%f','Delimiter',',');
    tline = fgets(fid);
    data = [data; C{1}'];
end
fclose(fid);
%end section taken from Andy's code
%textdata is the headings which contains the element names
%data is the numeric data

figure (fignum)
if strcmp(clearfig,'Yes')==1
clf
end

%extract the data wanted%
Nelements = size(textdata,2) - 5; %minus 5 due to XYZ, Unranged and radiusgyration%
Nclusters = size(data,1);
NnumerIon = zeros(Nclusters,1);
NdenomIon = zeros(Nclusters,1);
for i = 1:Nelements + 5
for ion = 1:size(numerIon,2)
    if strcmp(textdata{i},numerIon{ion}) == 1
        NnumerIon = NnumerIon + data(:,i); %number of numerator ions in each cluster%
    end
end
for ion = 1:size(denomIon,2)
    if strcmp(textdata{i},denomIon{ion}) == 1
        NdenomIon = NdenomIon + data(:,i); %number of denomiator ions in each cluster%
    end
end
end
%determine size of clusters%
if strcmp(xvalue, 'rg') ==1
    clustersize = data(:,Nelements + 5); %radius of gyration values%
elseif strcmp(xvalue, 'nions') == 1
    clustersize = NdenomIon; %using the number of denominator ions as cluster size
elseif strcmp(xvalue, 'ratio') == 1
    clustersize = NdenomIon + NnumerIon; %using the sum of ions in the ratio as cluster size%
elseif strcmp(xvalue, 'Z') == 1
    clustersize = data(:, 3);
elseif size(xvalue,1) == Nclusters
    disp('plotting alternative x values')
    clustersize = xvalue;
else 
    disp('format of xvalue incorrect, must be rg, nions, ratio or a column vector of alternate values')
    return
end

%form bins%
if isscalar(nbins) == 1
binends = logspace(log10(min(clustersize)),log10(max(clustersize)),nbins+1);%NB binends(1) is actually a bin start%
binends(nbins+1) = binends(nbins+1) + 0.0000000000001; %makes end of last bin bigger by one?%
else 
    binends = nbins;
    nbins = size(binends,2) - 1;
end

plotdata = zeros(nbins,4);
if strcmp(denomIon,'one')== 1 %incase I want to plot just the number of ions, not a fraction or ratio
    NdenomIon = ones(Nclusters,1);
end
fraction = NnumerIon./NdenomIon; %The compositon or ratio being plotted%

%4 columns 
%1. bin mid
%2. no. clusters in the bin
%3. av. of fraction 
%4. error in fraction
for i = 2:nbins+1
    binnedF = zeros(Nclusters,1);
    numinbin = 0;
    for j = 1:Nclusters
        if  (clustersize(j)) < binends(i) && (clustersize(j) >= binends(i-1))
            binnedF(j,1) = fraction(j);
            numinbin = numinbin+1;
        else
            binnedF(j,1) = NaN;
        end
    end
    plotdata(i-1,1) = (binends(i)-binends(i-1))/2 + binends(i-1);%bin midpiont%
    plotdata(i-1,2) = numinbin; %number of clusters in that size range%
    plotdata(i-1,3) = mean(binnedF,'omitnan');% av. of fraction%
    if strcmp(erbar,'stddev') == 1 %plots error bar of 1 std devition%
        plotdata(i-1,4) = std(binnedF,'omitnan');
    else
        plotdata(i-1,4) = std(binnedF,'omitnan')/sqrt(numinbin); %standard error in mean%%assuming the pionts are normally distributed etc.%
    end
    if plotdata(i-1,2) == 0
        plotdata(i-1,2) = NaN;
    end
    if plotdata(i-1,4) == 0
        plotdata(i-1,4) = NaN;
    end
end
averagecomp = [mean(fraction,'omitnan'),std(fraction,'omitnan')/sqrt(Nclusters)] ;
%above here is where all the outputed numbers are formed%

%scaling of the spots%
if strcmp(scalespots,'No') == 1
    scalespots = ones(nbins,1)*linearscale;
elseif strcmp(scalespots,'Yes') == 1 %scalespot size by number of clusters in that bin%
    scalespots = plotdata(:,2)*linearscale;
    if plotdata(1,2) == 1
        label1 = ' precipitate';%could build in functionality to decide if it's a ppt or cluster%
    else 
        label1 = ' precipitates';
    end
    if plotdata(nbins,2) == 1
        label2 = ' precipitate';
    else 
        label2 = ' precipitates';
    end
    hold on
    %labels is scaling spot by number of precipitates/cluster%
    text(plotdata(1,1),plotdata(1,3),{num2str(plotdata(1,2)) label1 '\downarrow'},'HorizontalAlignment','center','VerticalAlignment','bottom')
    text(plotdata(nbins,1),plotdata(nbins,3),{num2str(plotdata(nbins,2)) label2, '\downarrow'},'HorizontalAlignment','center','VerticalAlignment','bottom')
elseif size(scalespots,2) == nbins
    disp('correct length of scalespots')
else
    disp('incorrect length of scalespots')
    return
end

%plot data%
hold on
scatter(plotdata(:,1),plotdata(:,3),scalespots,Colour,'filled')  
    hold on
    if strcmp(erbar,'stder') == 1
        errorbar(plotdata(:,1),plotdata(:,3),plotdata(:,4),'k.')
    elseif strcmp(erbar,'stddev') == 1 %plots error bar of 1 std devition%
        errorbar(plotdata(:,1),plotdata(:,3),plotdata(:,4),'k.')
    elseif strcmp(erbar, 'No') ==1 
        disp('not plotting error bar')
    else 
        disp('error bar plttoing options, stder, stddev or No') 
    end
    xlim([min(clustersize),max(clustersize)])
    %ylim([0,max(plotdata(:,3))])
    if strcmp(xvalue,'rg') == 1
        label3 = 'cluster size, radius of gyration / nm';
    else
        label3 = 'cluster size';
    end
    xlabel(label3,'fontsize',12)
    label4 = '(';
    for i = 1:size(numerIon,2)
        if i ==1
            label4 = [label4 numerIon{i}];
        else
            label4 = [label4 '+' numerIon{i}];
        end
    end
    label4 = [label4 ')/('];
    for i = 1:size(denomIon,2)
        if i ==1
            label4 = [label4 denomIon{i}];
        else
            label4 = [label4 '+' denomIon{i}];
        end
    end
    label4 = [label4 ')'];
    ylabel(['Cluster composition ' label4],'fontsize',12)
    set(gca, 'XScale', 'log') %comment out this line to go back to linear scale%
    %would also need to change limits a little as well%
    if strcmp(plotmean,'Yes') == 1
    plot([min(clustersize),max(clustersize)],[mean(fraction,'omitnan'),mean(fraction,'omitnan')],'k')
    end
    legend(label4,'standard error',['mean ' label4],'Location','North')
    hold off
end

