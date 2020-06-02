function [fraction] = COI_clustercomp2(counts,header,allCOIID,COIID,numerIon,denomIon) 
%ecletic boogaloo
%old function takes a clusterstats file, this one works with counts and
%header matrices outputed by clustercomp_from_indxclrpos

%Inputs%
%counts - matrix of counts where each row is a cluster
%header - identity of species in each column
%allCOIID - list of all IDS in the counts table, needed as the row number
%does not nesecarily equal the ID, it  is check that it's the correct length
%COIID - IDS of the clusters that are to have the composition calucalted
%for. Can set as 'all'
%numerIon - numerator ion.
%denomIon - denominator ion/ions for which you use to calculate composition
%or ratio. Can be set as 'one'.

%Outputs%

%Comments
%averaging problem?

if strcmp(COIID,'all') == 1
    COIID = allCOIID;
end

%extract the rows of interest
[~,ia,~] = intersect(allCOIID,COIID); %find indices of COIID in whole list
COIcounts = counts(ia,:);

Nelements = length(header); 
Nclusters = size(COIcounts,1); %number of COI , not total number
NnumerIon = zeros(Nclusters,1);
NdenomIon = zeros(Nclusters,1);
for i = 1:Nelements
    for ion = 1:size(numerIon,2)
        if strcmp(header{i},numerIon{ion}) == 1
            NnumerIon = NnumerIon + COIcounts(:,i); %number of numerator ions in each cluster%
        end
    end
    for ion = 1:size(denomIon,2)
        if strcmp(header{i},denomIon{ion}) == 1
            NdenomIon = NdenomIon + COIcounts(:,i); %number of denomiator ions in each cluster%
        end
    end
end

if strcmp(denomIon,'one')== 1 %incase I want to plot just the number of ions, not a fraction or ratio
    NdenomIon = ones(Nclusters,1);
end
fraction = NnumerIon./NdenomIon; %The compositon or ratio being calcualted%