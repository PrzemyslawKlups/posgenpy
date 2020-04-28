function [x,y,z,m,nb]=readpos(fileName)
%[x,y,z,m,nb]=readpos(fileName)
% function to read .pos files
% useage [x,y,z,m,nb]=readpos('fileName')
% this returns, for each ion, the x,y,z postion
% m/z ratio and multiple hit identity
% Original E. Marquis, modified A. London Oct 2012
%
%   AtomProbeLab  Copyright (C) 2017  Andrew J. London\n
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.\n
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.\n
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

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