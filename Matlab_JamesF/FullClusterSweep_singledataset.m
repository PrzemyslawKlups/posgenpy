%does a full sweep of dmax, producing and saving text files with data for 
%each value of dmax, both for real data and one relabelled comparator.
%Then plots various figure to aid cluster parameter selection%

% posfile = 'E:\NetBeansProjects\PFIB_4_v33_dPA_TMA6_M1_R5111_09252\recons\recon-v01\default\R5111_09252-v01.pos';
% rangefile ='E:\NetBeansProjects\PFIB_4_v33_dPA_TMA6_M1_R5111_09252\recons\recon-v01\default\R5111_09252-v01.rrng';
% outputfolder ='C:\Users\James\Documents\Analysis\ref&v33_Jun19\ClusterSweep\9252_test';
% posoutputfolder = 'E:\NetBeansProjects\PFIB_4_v33_dPA_TMA6_M1_R5111_09252\recons\recon-v01\default';
% samplename = 'R5111_09252-v01';
% posfile = 'E:\NetBeansProjects\v33_dPA_1_R5111_08837\recons\recon-v02\R5111_08837-v02-clip-noCutip\R5111_08837-v02-clip-noCutip.pos';
% rangefile ='E:\NetBeansProjects\v33_dPA_1_R5111_08837\recons\recon-v02\R5111_08837-v02-clip-noCutip\R5111_08837-v02-clip-noCutip.rrng';
% outputfolder ='C:\Users\James\Documents\Analysis\ref&v33_Jun19\ClusterSweep\8837_test';
% posoutputfolder = 'E:\NetBeansProjects\v33_dPA_1_R5111_08837\recons\recon-v02\R5111_08837-v02-clip-noCutip';
% samplename = 'R5111_08837-v02-clip-noCutip';
posfile = 'E:\NetBeansProjects\PFIB_1_v33_dPA_M6_R5083_08340\recons\recon-v01\default\R5083_08340-v01.pos';
rangefile ='E:\NetBeansProjects\PFIB_1_v33_dPA_M6_R5083_08340\recons\recon-v01\default\R5083_08340-v01-alpha.rrng';
outputfolder ='C:\Users\James\Documents\Analysis\ref&v33_Jun19\ClusterSweep\8340';
posoutputfolder = 'E:\NetBeansProjects\PFIB_1_v33_dPA_M6_R5083_08340\recons\recon-v01\default';
samplename = 'R5083_08340-v01';
% posfile = 'E:\NetBeansProjects\v33_T7_2_R83_06144\recons\recon-v02\default\R83_06144-v02.pos';
% rangefile ='E:\NetBeansProjects\v33_T7_2_R83_06144\recons\recon-v02\default\R83_06144-v02-alpha.rrng';
% outputfolder ='C:\Users\James\Documents\Analysis\ref&v33_Jun19\ClusterSweep\6144';
% posoutputfolder = 'E:\NetBeansProjects\v33_T7_2_R83_06144\recons\recon-v02\default';
% samplename = 'R83_06144-v02';
dinitial = 0.1;
dstep = 0.05;
dfinal = 2;
inclusion = 4; %rember doing an inclusion would change dmax and Nmin values I select, especially Nmin%
ClusteredIons = ["Mg", "Si", "Cu"];
BulkIons = ["Al","V","AlH","AlH2"];
inputfolder = outputfolder;
%Nmin is a clip once clusteres are formed, Nmin selection only effects plotting
Nminmin = 5; %min value of Nmin to plot
Nminmax = 40; %max value of Nmin to plot
Nminstep = 1; %stepsize in Nmin%
%Clustersweeping step
%% Cluster Sweeping step - time consuming, commemt out to use data already saved
% ClusterSweep_incStats(posfile,rangefile,ClusteredIons,BulkIons,dinitial,dstep,dfinal,inclusion,outputfolder)
% ClusterSweep_incStatsRL(posfile,rangefile,ClusteredIons,BulkIons,dinitial,dstep,dfinal,inclusion,outputfolder)
%% Parameter Selection and Plotting

% this bit is for cluster parameter selection%
p = 0.95; %fractional threshold of "real" clusters required.
sizedist_Struct = form_sizedist_Struct(dinitial,dstep,dfinal,inputfolder,4); %extract saved data%
sizedistRL_Struct = form_sizedistRL_Struct(dinitial,dstep,dfinal,inputfolder,4);

Total_clusters_func_dmax_Nmin(dinitial,dstep,dfinal,Nminmin,Nminmax,sizedist_Struct);

%see select parameter function for full details on plots produced. Recommended
[Ncdifhighp_Struct,Nchighp_Struct,percent_real_Struct] = select_parameters(dinitial,dstep,dfinal,Nminmin,Nminmax,sizedist_Struct,sizedistRL_Struct,p);

dif_func_dmax_Nmin_3D %the higlights the section of the figure 3 were dmax and Nmin are to be excepted, also plots this region in seperate figures.
%3D plot of composition of cluster as a function of dmax and Nmin
%% Composition as a function of dmax and Nmin
%takes a little longer to run so comment out if not using
numerIon = {'Cu'};
denomIon = {'Mg','Si','Cu'};
% [compositions, matrixcompositions,compositionaca,NnumerIonsClstr,NdenomIonsClstr] =  Comp_clusters_func_dmax_Nmin(0.1,dstep,dfinal,inputfolder,Nminmin,Nminmax,Nminstep,numerIon,denomIon,ClusteredIons,posfile,rangefile);

%% Highlight single value of dmax and Nmin
dmax = 0.7;
Nmin = 20;

fname = ['dmax' num2str(dmax*100)]; %field name
figure(3)
hold on
scatter3(Nmin,dmax,sum(sizedist_Struct.(fname)(Nmin:size(sizedist_Struct.(fname),1),2)),'filled','g') %green spot%
hold off

