function [COIID_clipped] = clipCOI(COIID, inputparameter, cutoff, direction,varargin)
%COI - IDS of clusters, can set as 'all' but this will then require
%additional information of varargin set as indxClrPos, otherwise this input isn't needed
%inputparameter - the parameter you are using to clip the clusters, i.e.
%could be sizenatoms, allaxes(:,4,1), rg or some other thing
%cutoff - the boundary value
%direction - 'lessthan' or 'morethan'
%varargin - indxClrPos, path to the index cluster posfile or left blank

%NB - currenlty does less than or more than, not equal to
%NB - inputparameter must be same length as COIID, however I have been
%using COIID has a row (1xn, matrix) and all the parameters I outputs are
%columns (nx1) so this could get a little confusing. As they're both single
%line matrices though I can compare length and  can index them by one value so currently this is
%worked around, but may cause fututre confusion

%could add in a between here as well as more than lessthan, or rewrite it
%even more to have it like a histogram - although could maybe also achieve
%this with a between functionality and for loop where I vary the window. 

narginchk(4,5)
if strcmp (COIID,'all') == 1 
    [~, ~, ~, m, ~] = readpos(varargin{1,1});
    COIID = unique(m);
end 
if length(COIID) == length(inputparameter) 
else
    disp('number of IDS supplied does not match with corresponding number of inputparameters')
    return
end

COIID_clipped = 0;
for i = 1:length(COIID)
    if strcmp(direction,'lessthan') == 1
        if inputparameter(i) < cutoff
            COIID_clipped = [COIID_clipped, COIID(i)];
        end    
    elseif strcmp(direction,'morethan') == 1
        if inputparameter(i) > cutoff
            COIID_clipped = [COIID_clipped, COIID(i)];
        end  
    else
        disp('incorrect direction')
    end
end
COIID_clipped = COIID_clipped(2:size(COIID_clipped,2));
COIID_clipped = double(COIID_clipped);
end