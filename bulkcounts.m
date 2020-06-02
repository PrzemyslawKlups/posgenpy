function [bulkcount,header] = bulkcounts(posfile,rangefile) 
%inputs%
%pos file
%a posfile which contains the mass/charge ratios of all the points i.e. could be the whole dataset, or the clustered
%dataset, depending on what you want
%rangefile to go with the non-indexed pos file%

%outputs
%bulkcout - list of counts of each ranged species
%header - names of ranged species

%Comments
%did have this function doing a fraction caluclation, that code is still
%here but I think unessecary. would need to alter the function line at the
%top to output it again
% ,fraction]   ,numerIon,denomIon)

%test input values
% posfile = 'C:\Users\James\Documents\NetBeansProjects\AA6008_iT7_1_R5083_07423\recons\recon-v02\R5083_07423-v02-roi-no-Cu-tips\R5083_07423-v02-roi-no-Cu-tips.pos';
% rangefile = 'C:\Users\James\Documents\NetBeansProjects\AA6008_iT7_1_R5083_07423\recons\recon-v02\R5083_07423-v02-roi-no-Cu-tips\R5083_07423-v02-alpha.rrng';
% % numerIon = {'Mg'};
% denomIon = {'Mg', 'Si', 'Cu'};

DOMnode = xmlread(['bulkcounts.xml']);%loads xml file which is already in correct input format%
DOMnode.getElementsByTagName('posload').item(0).setAttribute('file',posfile); %Load indexed file%
DOMnode.getElementsByTagName('composition').item(0).setAttribute('rangefile',rangefile); %range file for mass to charge data%
xmlwrite('posgeninput_bulkcounts.xml',DOMnode);%write final xml fed to posgen%
!posgen.exe & ;
[status,cmdout] = dos('posgen posgeninput_bulkcounts.xml'); %excutes posgen%
%note by using cmdout the text will no longer appear in the command window,
%could be useful elsewhere also%
k = strfind(cmdout,'Counts');
data = cmdout(k+7:size(cmdout,2));
k = strfind(data,'Unranged');
data = data(1:k+22);
[bulkcount, na, errmsga] = sscanf(data,'%*s%d');
header = cell(size(bulkcount));
%str = string(data); %if I want output to be strings, dunno what I want yet
str = data;
for i = 1:size(bulkcount)
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
end