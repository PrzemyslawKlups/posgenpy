function [x,y,z,m,nb]=readpos(fileName)
%[x,y,z,m,nb]=readpos(nom_fichier)
% function to read .pos files
% useage [x,y,z,m,nb]=readpos('fileName')
% this returns, for each ion, the x,y,z postion
% m/z ratio and multiple hit identity
% Original E. Marquis, modified A. London Oct 2012
name = regexprep(fileName,'\.pos',''); % removes possible extention
fix='.pos';
name=strcat(name,fix); % adds file extention
[fid,msg] = fopen(name, 'r');
if fid==-1 % error checking
    disp('error');
    error([msg ' ' name]);
end
disp('Reading pos file...');

lflo=fread(fid, [4,inf], '*float32', 'b');
nb=size(lflo,2);
%flo=reshape(lflo,[4 nb]);
[x,y,z,m]=deal(lflo(1,:),lflo(2,:),lflo(3,:),lflo(4,:));
clear lflo
fclose(fid);
disp(strcat('OK, finished reading:',fileName,':',num2str(nb),' ions'));

end