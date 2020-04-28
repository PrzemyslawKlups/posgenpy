function [element_num, range_num, elements, ranges] = rangeReader (rangeFile)
% [element_num, range_num, elements, ranges] = rangeReader (rangeFile)
% reads a rrng file of the format:
% [Ions]
% Number=9
% Ion1=Fe
% Ion2=C
% Ion3=N
% ...
% Ion9=W
% [Ranges]
% Number=59
% Range1=27.9475 27.9913 Vol:0.01177 Fe:1 Color:FF00FF
% Range2=55.8950 55.9738 Vol:0.01177 Fe:1 Color:FF00FF
% ...
% Range24=69.9125 69.9738 Vol:0.04060 Fe:1 O:1 Color:FF0000
% Range59=23.9750 24.0363 Vol:0.01757 C:2 Color:660066
%
%
% Returns: number of elements and ranges, cell list of elements and cell
% containing the range in a readable format:
% 1          2        3   4      5            4+element_num
% rangeStart rangeEnd vol colour element1 ... elementN
%
% Where element1 is the number of elements in that ion in that range
% open and read rrng file
% WARNING has some hard coded working to reorder the range file (line55)
% Known bug: If the number of ions is wrong then it will fail to open the
% range file, or throw some other random error
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

verbose = 0; % display progress? 1/0
if exist(rangeFile,'file')
fid = fopen(rangeFile);
else
    error(strcat('Range file:',rangeFile,' not found'));
end
tline = fgets(fid); % this is the first line from the rrng file, [Ions]
if(strcmp(cellstr(tline),'[Ions]'))
    % note the cellstr removes the white space
    %keep going
    if verbose
        disp('Reading ions...');
    end
else
    error('Error reading [Ions] tag from range file');
    return;
end
% get the 2nd line
tline = fgets(fid);
% this is the number of elements
element_num = sscanf(tline,'Number=%f');
% Display number of elements
if verbose 
disp(strcat('Number of elements=',num2str(element_num)));
end
%make the elements cell array
elements = cell(element_num,1);
% for the next element_num lines, read the ions into the elements list
for e = 1:element_num
    tline = fgets(fid);
    S =  textscan(tline,'Ion%*u8=%s');
    elements(e,1) = S{1};
end

% check the order of the elements
elementIdeal = {'Fe' 'Cr' 'Y' 'Ti' 'W' 'Cu' 'Mo' 'Nb' 'Cd' 'Se' 'Ga' 'Ba' 'Ga' 'S' 'O' 'C' 'N' 'H'}';
[elements,~] = elementSort(elements,elementIdeal);

% so Ion1 is elements{1} = 'Fe' for example
% the next line should be [Ranges]
tline = fgets(fid);
if(strcmp(cellstr(tline),'[Ranges]'))
    % note the cellstr removes the white space
    %keep going
    if verbose 
    disp('Read elements');
    end
else
    disp('Error reading elements from range file');
    return;
end
% get the next line
tline = fgets(fid);
% this is the number of ranges
range_num = sscanf(tline,'Number=%f');
extraCols = 4;
ranges = cell(range_num,element_num+extraCols); % make an array to store the range details

for r = 1:range_num
    tline = fgets(fid);
    C = textscan(tline,'Range%*u8=%f %f');
    ranges{r,1} = C{1};
    ranges{r,2} = C{2};
    % find the elements present in ranged ion
    for e = 1:element_num
        % this adds the number of each element present to the 'ranges' table
        % it scans for each element, e, in the line from the range file
        % then it adds on the right number of characters to move to the
        % position of the number of that element present in the ranged ion
        % The strcat...':' is so that it finds 'Co:#' not 'Color ect'
        % disp(strfind(tline,strcat(elements{e},':')));
        numOfE = str2double(tline(strfind(tline,strcat(elements{e},':'))+length(elements{e})+1));
        if(~isnan(numOfE))
            ranges{r,extraCols+e} = numOfE;
        else
            ranges{r,extraCols+e} = 0;
        end
    end
    % find ion volume information
    vol = textscan(tline,'%*s %*f Vol:%f');
    ranges{r,3} = vol{1};
    % find colour information
    color = textscan(tline(strfind(tline,'Color:'):end),'Color:%s');
    ranges(r,4) = color{1};
end
fclose(fid);
end