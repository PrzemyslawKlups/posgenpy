clear
%% Input data, posfiles and any poles if using 
%6195
indxClrPos = 'E:\NetBeansProjects\v33_T7_4_R83_06195\recons\recon-v01\default\R83_06195-v01.cluster.indexed.pos';
masschargedata = 'E:\NetBeansProjects\v33_T7_4_R83_06195\recons\recon-v01\default\R83_06195-v01.cluster.pos';
rangefile = 'E:\NetBeansProjects\v33_T7_4_R83_06195\recons\recon-v01\default\R83_06195-v01.rrng';
poles_detector_coord = 'C:\Users\James\Documents\MATLAB\Crystalography\poles6195.csv'; %matrix of hkl values and detector postion%
%% Example of fitting all axes
%ploting all clusters morphology and stereoplot of orienation of longest
%principle axes

[allaxes, sizenatoms, clustercenter] = Fit_principle_axes(indxClrPos,1,'No');
plotdataspatter = morphologygraph(allaxes,sizenatoms,2,'r',0.1,'detected atoms','Yes');
plotdatastereo = stereographplot(allaxes,sizenatoms,3,'Yes','No','long','default',0.01);

%Keys%
%plotdata_morph =  morphologygraph(allaxes,scalespots,fignum,colour,linearscale,label,clearfig)
%plotdata_stereo = stereographplot(allaxes,scalespots,fignum,clearfig,rotated,WhichAxes,Colours,linearscale)

%% Transforming axes to crystal space

ICF = 'Matlab';
[polestipframe,R,ICFout] = transformaxes(poles_detector_coord,ICF,42); %I assume Xand Y on detector are in same direction as X and Y in reconstruction, this may not be the case%
%virtual fight path of XR is 42. XS it's ____
[rotatedallaxes,integerdirections] = rotateaxes(allaxes,R);
plotdata_stereo2 = stereographplot(rotatedallaxes,sizenatoms,4,'Yes','Yes','long','default',0.01);

%% Adding XYZ to crystal space plot or <001> to reconstruction space plot
tipXYZincrys = zeros(1,4,3);
tipXYZincrys(1,1:3,:) = R; %direction are just the column vectors in R (I*R = R and identiy matrix is the X, Y and Z axes%
tipXYZincrys(1,4,:) = 1;
crysoooneinXYZ = zeros(1,4,3); %crystal <001> direction in detector reconstruction coordinates%
crysoooneinXYZ(1,1:3,:) = inv(R); %direction are just the column vectors in R (I*R = R and identiy matrix is the X, Y and Z axes%
crysoooneinXYZ(1,4,:) = 1;
XYZstereo = stereographplot(tipXYZincrys,tipXYZincrys(:,4,1),4,'No','Yes','longmidshort',{'m','c','k'},50);
ooonestereo = stereographplot(crysoooneinXYZ,crysoooneinXYZ(:,4,1),3,'No','No','longmidshort',{'c','c','c'},50);

%% Higlightling one cluster

COIID = [1];
[allaxes_COI, sizenatoms_COI, clustercenter_COI,P] = Fit_principle_axes_subset(indxClrPos,5,COIID,'Yes');
plotdata_morph_COI =  morphologygraph(allaxes_COI,sizenatoms_COI,2,'g',0.1,'No','No');
plotdatastereo_COI = stereographplot(allaxes_COI,sizenatoms_COI,3,'No','No','longmidshort',{'g',[0,0.5,0],[0,0.25,0]},0.01);

%% Histogram of axis lengths

figure(6)
clf
hold on
title('histogram of axis lengths')
hist(allaxes(:,4,:))
hold off

%% Scale spots by radius of gyration or composition

[rg,allCOIID] = rggenerator(indxClrPos,'all');
plotdata_morph_COI =  morphologygraph(allaxes,rg,7,'r',10,'rg','Yes');

% [counts,header] = clustercomp_from_indxclrpos(indxClrPos,masschargedata,rangefile,'all');
% [fraction] = COI_clustercomp2(counts,header,allCOIID,'all',{'Cu'},{'Cu','Mg','Si'});
% plotdata_morph_COI =  morphologygraph(allaxes,fraction,8,'r',10,'rg','Yes');

%% Clip out smaller ppts

[COIID_clipped] = clipCOI('all', sizenatoms,1000,'morethan',indxClrPos);
[allaxes_COI_clipped, sizenatoms_COI_clipped, clustercenter_COI_clipped,P] = Fit_principle_axes_subset(indxClrPos,8,COIID_clipped,'No');
plotdata_morph_COI_clipped =  morphologygraph(allaxes_COI_clipped,sizenatoms_COI_clipped,2,'b',0.1,'No','No');

