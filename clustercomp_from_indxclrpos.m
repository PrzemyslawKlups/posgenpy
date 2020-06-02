function [counts,header] = clustercomp_from_indxclrpos(indxClrPos,masschargedata,rangefile,COIID) 
%inputs%
%indexed cluster pos file (probably after a splitting of clusters has been
%done)
%a posfile which contains the mass/charge ratios of all the points in the
%indexed cluster file, i.e. could be the whole dataset, or the clustered
%dataset, clustered dataset will be quicker.
%rangefile to go with the non-indexed pos file%
%COIID - IDS of clusters you want the composition for, can set as 'all'

%outputs
% counts- each row is the counts for one clsuter, each column is one ranged
% species
%header - names of the ranged species, in same order as columns of counts

%Comments
%did have this function doing a fraction caluclation, that code is still
%here but I think unessecary. would need to alter the function line at the
%top to output it again
% ,fraction]   ,numerIon,denomIon)

%test input values
% indxClrPos = 'C:\Users\James\Documents\NetBeansProjects\AA6008_iT7_1_R5083_07423\recons\recon-v02\R5083_07423-v02-roi-no-Cu-tips\R5083_07423-v02-roi-no-Cu-tips-alpha-dmax-0.9-Nmin-17-split.cluster.indexed.pos';
% masschargedata = 'C:\Users\James\Documents\NetBeansProjects\AA6008_iT7_1_R5083_07423\recons\recon-v02\R5083_07423-v02-roi-no-Cu-tips\R5083_07423-v02-roi-no-Cu-tips-alpha-dmax-0.9-Nmin-17-split.cluster.pos';
% rangefile = 'C:\Users\James\Documents\NetBeansProjects\AA6008_iT7_1_R5083_07423\recons\recon-v02\R5083_07423-v02-roi-no-Cu-tips\R5083_07423-v02-alpha.rrng';
% COIID = 'all';
% numerIon = {'Mg'};
% denomIon = {'Mg', 'Si', 'Cu'};

%read pos indx clr to get IDs, or have manual option for IDS
% loop, write range file, pass to posgen, extract composition, add to final
% data

%write range file, just need to write a text file of correct format and
%make the suffix .rrng . I think. 

if strcmp (COIID,'all') == 1
    [~, ~, ~, m, ~] = readpos(indxClrPos);
    IDS = unique(m);
else
    IDS = COIID;
end  

counts = 1; %don't yet know the number of ranges so waiting till end of first loop to preallocate%
for row = 1:size(IDS,2)  %there might be missing values in a list of IDS, 
    %the output matrix "counts" will not have gaps and the clusters will be 
    %the in the order specified by IDS (which is not nessacrily in numeric
    %order
ID = IDS(row); %solves the problem described just above%
%write rangefile%
indxrangefile = 'singlecluster_index.rrng';
lowerbound = ID - 0.1000;
upperbound = ID + 0.1000;
fileID = fopen(indxrangefile,'w');
fprintf(fileID,'%6s\n','[Ions]');
fprintf(fileID,'%7s\n','Number=1'); 
fprintf(fileID,'%13s\n','Ion1=Cluster1');
fprintf(fileID,'%7s\n','[Ranges]');
fprintf(fileID,'%7s\n','Number=1');
fprintf(fileID,'%7s%6.4f %6.4f %38s','Range1=',lowerbound,upperbound,'Vol:0.00000 Name:Cluster1 Color:0000FF'); 
fclose(fileID);


%then can use code below to make the new pos files%
DOMnode = xmlread(['clustercomp_from_indxclrpos.xml']);%loads xml file which is already in correct input format%
DOMnode.getElementsByTagName('posload').item(0).setAttribute('file',indxClrPos); %Load indexed file%
DOMnode.getElementsByTagName('range').item(0).setAttribute('file',indxrangefile);%range file for index, giving just one cluster%
DOMnode.getElementsByTagName('posload').item(1).setAttribute('file',masschargedata); %Could maybe speed up the intersection by interesecting with the clustered file not the whole dataset 
DOMnode.getElementsByTagName('composition').item(0).setAttribute('rangefile',rangefile); %range file for mass to charge data%
xmlwrite('posgeninput_ClrCompfromIndxClrPos.xml',DOMnode);%write final xml fed to posgen%
!posgen.exe & ;
[status,cmdout] = dos('posgen posgeninput_ClrCompfromIndxClrPos.xml'); %excutes posgen%
%note by using cmdout the text will no longer appear in the command window,
%could be useful elsewhere also%
k = strfind(cmdout,'Counts');
data = cmdout(k+7:size(cmdout,2));
k = strfind(data,'Unranged');
data = data(1:k+22);
[snglclrcounts, na, errmsga] = sscanf(data,'%*s%d');
if counts == 1
    counts = zeros(size(IDS,2),size(snglclrcounts,1));
end
counts(row,:) = snglclrcounts;
end
    %just want it to get the headers once so doing this outside the loop, to
    %the last cluster
header = cell(size(snglclrcounts));
%str = string(data); %if I want output to be strings, dunno what I want yet
str = data;
for i = 1:size(snglclrcounts)
    [token, remain] = strtok(str);
    header{i,1} = token;
    [token, remain] = strtok(remain);
    str = remain;
end

%output a single composition or ratio%
% Nelements = size(header,2); 
% Nclusters = size(counts,1);
% NnumerIon = zeros(Nclusters,1);
% NdenomIon = zeros(Nclusters,1);
% for i = 1:Nelements + 5
% for ion = 1:size(numerIon,2)
%     if strcmp(header{i},numerIon{ion}) == 1
%         NnumerIon = NnumerIon + counts(:,i); %number of numerator ions in each cluster%
%     end
% end
% for ion = 1:size(denomIon,2)
%     if strcmp(header{i},denomIon{ion}) == 1
%         NdenomIon = NdenomIon + counts(:,i); %number of denomiator ions in each cluster%
%     end
% end
% end
% fraction = NnumerIon./NdenomIon; %The compositon or ratio being plotted%