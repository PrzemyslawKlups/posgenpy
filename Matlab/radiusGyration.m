function [rg,rx,ry,rz] = radiusGyration(x,y,z,m)
%[rg,rx,ry,rz] = radiusGyration(x,y,z,m)
% This function calculates the radius of gyration of points with masses m
% A. London August 2013
xm = sum(x.*m)/sum(m);
ym = sum(y.*m)/sum(m);
zm = sum(z.*m)/sum(m);
rx = sqrt((sum((x-xm).^2))/length(x));
ry = sqrt((sum((y-ym).^2))/length(y));
rz = sqrt((sum((z-zm).^2))/length(z));
rg = sqrt(rx*rx+ry*ry+rz*rz);
end