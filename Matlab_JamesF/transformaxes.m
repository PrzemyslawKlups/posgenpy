function [polestipframe,R,ICFout] = transformaxes(poles_detector_coord,ICF,flightpath)
%Transformaxes Takes xy coordinates and hkl values of poles in detector
%spcae, creates, xyz vectors of the poles in reconstruction space and then
%finds rotation matrix to tranform between reconstruction space and
%crystal space%
%   Detailed explanation goes here
%flightpath = 42; %hardcoded LEAP flight path, needed to caluclate angles. NB this is the virtual flight path for reflectron fitted LEAP5000HR%
%flightpath = 100; %supposedly for LEAP 5000 XS
%poles_detector_coord - file path to a csv file containing the xy position
%of the pole on the detector and an hkl assignation. Col 1 h, Col 2 k Col 3 l, Col
%4 x Col 5 y
%ICF - from IVAS 
%or my own calculations 
%or ICF = 'Matlab' to get this script to calcualte it, same method as my
%spreadsheet, average ICF of all poles combinations%
%R - rotation matrix, found by sinngle value decomposition solving of
%Wahba's problem (see below). 
%direction in crystall coordinates = R*direction in reconstuction coordinates

%code%
%load information%
detector = csvread(poles_detector_coord);
if size(detector,2) == 5
else
    disp('Incorrect number of columns in csv file');
   return
end

hkli = [detector(:,1), detector(:,2), detector(:,3)];
if strcmp(ICF,'Matlab') == 1
    %Calculate ICF%%Labbook 1 15/1/18 & 4/2/19%
    xi = detector(:,4);
    yi = detector(:,5); %the Geordie factor%
    ri = (xi.^2 + yi.^2).^0.5;%distance of pole from detector centre%
    phii = atan2(yi,xi); %azimuthal angle, anticlockwise from positive x-axis%
    thetai = atan(ri./flightpath); %angle between centre axis of tip and detector and path of tip to spot on detector%
    cosinethetacrysij = (hkli*hkli')./(vecnorm(hkli,2,2)*vecnorm(hkli,2,2)');%set of angles between two pols in the crystal%
    cosinethetaobsij = zeros(size(xi,1)); %set of angles between poles on the detector%
    for i = 1: size(xi,1)
        for j = 1:size(xi,1)
            if i == j
                cosinethetaobsij(i,j) = NaN;
                cosinethetacrysij(i,j) = NaN;
            else
                cosinethetaobsij(i,j) = cos(thetai(i))*cos(thetai(j)) + sin(thetai(i))*sin(thetai(j))*cos(phii(i)-phii(j));
            end
        end
    end
    ICFij = acos(cosinethetacrysij)./acos(cosinethetaobsij);
    ICF = sum(sum(ICFij,1,'omitnan'),2)/(size(ICFij,1)^2 - size(ICFij,1)) ;
end
ICFout = ICF;

%create vector in reconstruction coordiates using one value of ICF%
D = (detector(:,4).^2 + detector(:,5).^2).^0.5;%distance of pole from detector centre%
thetacrys = atan(D./flightpath).*ICF;  %polar angle of the crystallographic direction from tip centre axis, corrected by image compression%
phi = atan2(detector(:,5),detector(:,4)); %azimuthal angle anticlockwise from positive x%
polestipframe = [sin(thetacrys).*cos(phi),sin(thetacrys).*sin(phi),cos(thetacrys)]; 
%did wonderf if z should be a negative here as IVAS plots the tip of the reconstruction as zero, and doesn't change X and Y%
%but susqunet testing showed this gave a reflection aswell as rotation
%between poles in tip space and cystal space.

%Solve Wahba's problem with a single value decomposition%
%as described on wikipedia%
w = (hkli./vecnorm(hkli,2,2))';%normalise the vector%%no need to normalise crystal axes vectors, as vector magnitude doesn't make a difference?%
B = zeros(3,3); 
for i = 1:size(polestipframe,1)
a = w(:,i)*polestipframe(i,:); 
B = B + a;
end
[U,S,V]= svd(B); %single value decomposition% %function svds has more options, svd is not a one solution problem. svds fun
M = diag([1, 1, det(U)*det(V)]);
R = U*M*V';
end