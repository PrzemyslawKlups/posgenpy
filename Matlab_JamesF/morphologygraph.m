function [plotdata] =  morphologygraph(allaxes,scalespots,fignum,colour,linearscale,label,clearfig)
%Plots aspect ratio vs oblatness

%Outputs%
%plotdata = [aspectratio, oblateness, scalespots];

%Inputs%
%allaxes - n by 4 by 3 matrix of sxes and lengths, allready sorted by
%length(done automatically by fit axes function)
%scalespots - set to 'No' to not use, else it's a n by 1 matrix of some
%factor you want to scale the spot size by, where the order corresponds to
%the order in allaxes. eg allaxes(:,4,1) scale size by length of longest
%principle axis, or use sizenatoms from fit_principle_axes function to
%scale by number of atoms in cluster.
%fignum - number for figure, so figures aren't overwritten%
%colour - clour of marker as used by matlab scattter plot%
%linearscale - can scale he size of the cirlce so they better fit the
%graph, while maintaining the same scalling between different sets of data%
%value of 0.1 works well for 70000 ion ppt%
%label - adds in a label for scalespots function, lablels smallest and
%largest spots with the numerical value and the text specified here. Or
%set as 'No' to not label.
%clearfig - default yes so you're not writing onto a figure which already
%has stuff on, but can turn off to plot multple data sets on same graph%
figure (fignum)
if strcmp(clearfig,'Yes')==1
clf
end
hold on 
title('Precipitate Morphology L_1>L_2>L_3')
numberclusters = size(allaxes,1);
lengths = zeros(numberclusters,3);
for i =1:3
lengths(:,i) = allaxes(:,4,i);
end
aspectratio = lengths(:,2)./lengths(:,1);
oblateness = lengths (:,3)./lengths(:,2);
marceauplot = [aspectratio, oblateness, scalespots];
marceauplot = sortrows(marceauplot,3,'descend'); %sorting clusters by size%
scatter(marceauplot(:,2),marceauplot(:,1),marceauplot(:,3)*linearscale,'ok','MarkerFaceColor',colour)
xlabel('Oblateness L_3/L_2')
ylabel('Aspect ratio L_2/L_1')
text(0.9,0.95,'Sphere')
text(0.9,0.05,'Rod')
text(0.05,0.05,'Lath')
text(0.05,0.95,'Disc')
smallestppt = num2str(marceauplot(numberclusters,3));
largestppt = num2str(marceauplot(1,3));
if strcmp(label,'No') == 1
    disp('not labelling smallest and largest spots')
else  
text(marceauplot(numberclusters,2),marceauplot(numberclusters,1),{'\uparrow' smallestppt label }...
    ,'HorizontalAlignment','center','VerticalAlignment','top','Color','blue');
text(marceauplot(1,2),marceauplot(1,1),{'\uparrow' largestppt label }...
    ,'HorizontalAlignment','center','VerticalAlignment','top','Color','blue');
end
% plot([-2,2],[0.5,0.5],'-k') %no longer need this really but left just in
% case, noew using grid line instead of plotting in a black line%
% plot([0.5,0.5],[-2,2],'-k')
xlim([0,1]);
ylim([0,1]);
ax = gca;
ax.Position = [0.1,0.13, 0.8, 0.8]; %left bottom width height%
%leaving a space for  legen on the right%
ax.XTick = [0,0.5,1];
ax.XGrid = 'on';
ax.YTick = [0,0.5,1];
ax.YGrid = 'on';
ax.GridAlpha = 1;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.PlotBoxAspectRatio = [1,1,1];
hold off
plotdata = marceauplot; %I originally called it this as it was in a Marcaueau paper where I found this graph originally% 
end

