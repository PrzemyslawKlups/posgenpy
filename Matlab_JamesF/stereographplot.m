function plotdata = stereographplot(allaxes,scalespots,fignum,clearfig,rotated,WhichAxes,Colours,linearscale)
%stereographplot Plots a stereographic projection of axis directions
%
%Inputs%
%allaxes - n by 4 by 3 matrix of axes and lengths, allready sorted by
%length(done automatically by fit axes function)
%scalespots - set to 'No' to not use, else it's a n by 1 matrix of some
%factor you want to scale the spot size by, where the order corresponds to
%the order in allaxes. eg allaxes(:,4,1) scale size by length of longest
%principle axis, or use sizenatoms from fit_principle_axes function to
%scale by number of atoms in cluster.
%fignum - number figure to make, so you don't overwrite figures
%clearfig - whether or not to clear figure, turn of to allow multiple
%plotting on one graph.
%rotated - if yes changes axes lables from XYZ to [001] etc
%WhichAxes -string containing long mid and/or short as too which set of
%axes you want to plot. 
%Colours - 1 by 3 cell array containg the colours in which to plot the dots
%,or 'default', which just uses {'r','g','b'}%
%linearscale - can scale he size of the cirlce so they better fit the
%graph, while maintaining the same scalling between different sets of data%
%value of 0.1 works well for 70000 ion ppt. Has no effect if sizenatoms not used%

%convert to stereographic projections
%xs, x sterographic, ys, y stereograhic%
%from wolfram 
%xs = x/(1+z)
%ys = y(1+z)
%for negative values of Z though need to project from the other pole so it
%1 + absolute value of z
% L1 > L2 > L3
%convert all vectors so z value is positive - NB the three axis vector are
%still orthogonal they just will no longer be a righthanded set? if they
%ever were to begin with%
for i= 1:size(allaxes,1)
    for j = 1:3 
    if allaxes(i,3,j) < 0
        allaxes(i,1:3,j) = allaxes(i,1:3,j)*-1;
    end
    end
end
xs1 = allaxes(:,1,1)./(1+(allaxes(:,3,1)));
ys1 = allaxes(:,2,1)./(1+(allaxes(:,3,1)));
xs2 = allaxes(:,1,2)./(1+(allaxes(:,3,2)));
ys2 = allaxes(:,2,2)./(1+(allaxes(:,3,2)));
xs3 = allaxes(:,1,3)./(1+(allaxes(:,3,3)));
ys3 = allaxes(:,2,3)./(1+(allaxes(:,3,3)));
figure(fignum)
if strcmp(clearfig,'Yes')==1
clf
end
if strcmp(scalespots,'No') == 1
    scalespots = ones(size(allaxes,1),1);
    linearscale = 10;
end
if strcmp(Colours,'default')== 1
    Colours = {'r','g','b'};
elseif (iscell(Colours) == 1) && (size(Colours,2) == 3) && (size(Colours,1) == 1)
else
    disp('Incorrect format of Colours input, should be 1x3 cell array, or "default"')
    return
end
hold on
title({'Sterographic projection of axis vectors';' '})
if strfind(WhichAxes,'long') > 0
    %colour1 = ['.' Colours{1,1}];
    %plot(xs1,ys1,colour1)
    colour1 = Colours{1,1};
    l1 = [xs1,ys1, scalespots];
    l1 = sortrows(l1,3,'descend'); %sorting clusters by whatever factor used so smallest spots aren't hidden at the back%
    scatter(l1(:,1),l1(:,2),l1(:,3)*linearscale,'ok','MarkerFaceColor',colour1)
else
    l1 = zeros(size(allaxes,1),3);
end
if strfind(WhichAxes,'mid') > 0 
    colour2 = Colours{1,2};
    l2 = [xs2,ys2, scalespots];
    l2 = sortrows(l2,3,'descend'); %sorting clusters by whatever factor used so smallest spots aren't hidden at the back%
    scatter(l2(:,1),l2(:,2),l2(:,3)*linearscale,'ok','MarkerFaceColor',colour2)
else
    l2 = zeros(size(allaxes,1),3);
end
if strfind(WhichAxes,'short') > 0
    colour3 = Colours{1,3};
    l3 = [xs3,ys3, scalespots];
    l3 = sortrows(l3,3,'descend'); %sorting clusters by whatever factor used so smallest spots aren't hidden at the back%
    scatter(l3(:,1),l3(:,2),l3(:,3)*linearscale,'ok','MarkerFaceColor',colour3)
else
    l3 = zeros(size(allaxes,1),3);
end
theta = linspace(0,2*pi);
xcirc = cos(theta);
ycirc = sin(theta);
plot(xcirc, ycirc, '-k') %bounding cirlce%
if strcmp(rotated, 'Yes') == 1
    tx = text(1.1,0,'[ 1 0 0]');
    ty = text(-0.1,1.05,'[0 1 0]');
    plot([0,0],[-1,1],'-k')
    plot([-1,1],[0,0],'-k')
elseif strcmp(rotated, 'No') == 1
    tx = text(1.1,0,'X');
    ty = text(0,1.05,'Y');
else
    disp('incorrect form of input "rotated", should be "Yes" or "No"')
    return
end
xlabel(' L_1 > L_2 > L_3')
xlim([-1,1])
ylim([-1,1])
axis equal
hold off
plotdata = [l1,l2,l3];
end

