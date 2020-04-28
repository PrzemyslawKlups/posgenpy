function [rg,rx,ry,rz] = radiusGyrationMassless(x,y,z,m)
%[rg,rx,ry,rz] = radiusGyrationMassless(x,y,z,m)
% This function calculates the radius of gyration of points without masses m
% A. London August 2013
xm = mean(x);
ym = mean(y);
zm = mean(z);
rx = sqrt((sum((x-xm).^2))/length(x));
ry = sqrt((sum((y-ym).^2))/length(y));
rz = sqrt((sum((z-zm).^2))/length(z));
rg = sqrt(rx*rx+ry*ry+rz*rz);
end