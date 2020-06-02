function [COIalldata,header,COIcomp,COIcompav,COIcompstder,COItotalcomp] = COI_clustercomp(clusterstatsfile,COIID,numerIon,denomIon)
%Inputs%
%clusterstatsfile - clusterstats file ouputted by posgen- NEEDS TO BE
%ALREADY FILTERED BY NMIN ALREADY
%COIID - cluster/s of interest ID. 'all' gives all the data
%numerIon - numerator ion.
%denomIon - denominator ion/ions for which you use to calculate composition
%or ratio

%Outputs
%COIalldata - complete rows from the cluster stats file for the clusters of
%interest inclusing XYZ and rg always as counts
%header - header for COIalldata
%COIcomp - composition/ratio defined by numerion and denomion
%COIcompav - average of the compositions
%COIcompstder - standard error in the composition
%COItotalcomp - will give the composition of the
%cluster summeed together.


% clusterstatsfile = 'C:\Users\James\Documents\NetBeansProjects\v33_iPA_170_30m_1_R83_06453\recons\recon-v02\default\R83_06453-v02-clusterstats-dmax-0.7-Nmin-15.txt';
% COIID = [7,97];
% numerIon = {'Mg'};
% denomIon = {'Mg', 'Si', 'Cu'};

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
header = C{1}';
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

if strcmp(COIID,'all')==1
    COIID = 1:1:size(data,1);
end
%extract the COI wanted%
COIalldata = data(COIID,:);
%extract the data wanted%
Nelements = size(header,2) - 5; %minus 5 due to XYZ, Unranged and radiusgyration%
Nclusters = size(COIalldata,1);
NnumerIon = zeros(Nclusters,1);
NdenomIon = zeros(Nclusters,1);
for i = 1:Nelements + 5
for ion = 1:size(numerIon,2)
    if strcmp(header{i},numerIon{ion}) == 1
        NnumerIon = NnumerIon + COIalldata(:,i); %number of numerator ions in each cluster%
    end
end
for ion = 1:size(denomIon,2)
    if strcmp(header{i},denomIon{ion}) == 1
        NdenomIon = NdenomIon + COIalldata(:,i); %number of denominator ions in each cluster%
    end
end
end
COIcomp = NnumerIon./NdenomIon; %compostion for each COI%
COIcompav = mean(COIcomp,1); %average of these compositions%
COIcompstd = std(COIcomp,1); %standard deviation
COIcompstder =  COIcompstd/sqrt(Nclusters); %standard error%
COItotalcomp = sum(NnumerIon,1)/sum(NdenomIon,1); %Sum togehter cluster before calculating composition% 