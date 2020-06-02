function [compositions, matrixcompositions,compositionaca,NnumerIonsClstr,NdenomIonsClstr] = Comp_clusters_func_dmax_Nmin(dinitial,dstep,dfinal,inputfolder,Nminmin,Nminmax,Nminstep,numerIon,denomIon,ClusteredIons,posfile,rangefile)

%Composition of clusters as a function of dmax max and Nmin
%takes a bunch of clusterstats files produced via a cluster sweep using
%posgen and from that creates the cluster composition for a numerator and
%denominator elements.
% dintial, dstep and dfinal are dmax values that have been swept
% inputfolder is the folder in which the clusterstats files can be found
% Nminmin, smallest value of Nmin you wish to plot
% Nminmax largest value of Nmin which you wish to plot
%numerIon is the numerator ion, ie the ion your interested in. NB can be
%multiple Ions.
% denomIon is the other Ions in the cluster you want to compare the
% compoistion of numerIon to.
%ClusteredIons is the ions used to define the cluster in the original
%posgen sweep. This is needed to know which ions Nmin applies to.
%posfile is the whole dataset posfile and corresponding rangefile. They 
%are the needed for calulating matrix concentrations. 

[allbulk,header] = bulkcounts(posfile,rangefile);
bulkcountSOI = zeros(1,2); %Species of interest%%Will be total bulk numerator and denominator counts%
for j = 1:length(header)
    for i = 1:length(denomIon)
        if strcmp(denomIon{i},header{j})
            bulkcountSOI(1,2) = bulkcountSOI(1,2) + allbulk(j);
        end
    end
    for i = 1:length(numerIon)
        if strcmp(numerIon{i},header{j})
            bulkcountSOI(1,1) = bulkcountSOI(1,1) + allbulk(j);
        end
    end
end

T = readtable([inputfolder '\clusterstats-dmax-' num2str(dinitial) '.txt']);
header = T.Properties.VariableNames;
numerIoncols = zeros(1,size(numerIon,2));
for col = 1:size(header,2)
    for i = 1:size(numerIon,2)
        if strcmp(numerIon(i),header(col))
        numerIoncols(i) =  col;%columns which contains the numerIons%
        end
    end
end
denomIoncols = zeros(1,size(denomIon,2));
for col = 1:size(header,2)
    for i = 1:size(denomIon,2)
        if strcmp(denomIon(i),header(col))
        denomIoncols(i) =  col;%columns which contains the denomIons%
        end
    end
end
ClusteredIoncols = zeros(1,size(ClusteredIons,2));
for col = 1:size(header,2)
    for i = 1:size(ClusteredIons,2)
        if strcmp(ClusteredIons(i),header(col))
        ClusteredIoncols(i) =  col;%columns which contains the ClusteredIons%
        end
    end
end
numbIons = size(header,2) - 5;%Number of ion species that were in the range file supplied to posgen%
compositions = zeros((Nminmax-Nminmin)/Nminstep,round((dfinal-dinitial)/dstep));
compositionaca = zeros((Nminmax-Nminmin)/Nminstep,round((dfinal-dinitial)/dstep));
NnumerIonsClstr = zeros((Nminmax-Nminmin)/Nminstep,round((dfinal-dinitial)/dstep));
NdenomIonsClstr = zeros((Nminmax-Nminmin)/Nminstep,round((dfinal-dinitial)/dstep));
matrixcompositions = zeros((Nminmax-Nminmin)/Nminstep,round((dfinal-dinitial)/dstep));
for i = dinitial:dstep:dfinal %sweep in dmax, ie open each text file%
    clusterstats = dlmread([inputfolder '\clusterstats-dmax-' num2str(i) '.txt'],'',1,0);
    for Nmin = Nminmin:Nminstep:Nminmax
        NaboveNmin = 0;
        sumclustercomp = 0;
        NnumerIonsClustered = 0;
        NdenomIonsClustered = 0;
        for clusternum = 1:1:size(clusterstats,1)
           Nsoluteatoms = sum(clusterstats(clusternum,ClusteredIoncols));
           if Nsoluteatoms > Nmin
               clustercomp = sum(clusterstats(clusternum,numerIoncols))/sum(clusterstats(clusternum,denomIoncols));
               NaboveNmin = NaboveNmin + 1;
               sumclustercomp = sumclustercomp + clustercomp;
               NnumerIonsClustered = NnumerIonsClustered + sum(clusterstats(clusternum,numerIoncols));
               NdenomIonsClustered = NdenomIonsClustered + sum(clusterstats(clusternum,denomIoncols));
           end
        end
        compositions((Nmin-Nminmin)/Nminstep+1,round((i-dinitial)/dstep)+1) = sumclustercomp/NaboveNmin;%average composition of clusters, not the compostion of all the clusters ions%
        compositionaca((Nmin-Nminmin)/Nminstep+1,round((i-dinitial)/dstep)+1) = NnumerIonsClustered/NdenomIonsClustered; %overall composition of all the clustered ions%
        NnumerIonsClstr((Nmin-Nminmin)/Nminstep+1,round((i-dinitial)/dstep)+1)= NnumerIonsClustered; 
        NdenomIonsClstr((Nmin-Nminmin)/Nminstep+1,round((i-dinitial)/dstep)+1) = NdenomIonsClustered; 
        matrixcompositions((Nmin-Nminmin)/Nminstep+1,round((i-dinitial)/dstep)+1) = (bulkcountSOI(1,1) - NnumerIonsClustered)/(bulkcountSOI(1,2) - NdenomIonsClustered);
    end
end
[x,y]= meshgrid(Nminmin:Nminstep:Nminmax,dinitial:dstep:dfinal); %mesh of possible dmax and Nmin values values%
figure(7)
clf
mesh(x,y,transpose(compositions));
xlabel('N_m_i_n','fontsize',12)
ylabel('d_m_a_x','fontsize',12)
zlabel(['Cluster composition' numerIon{1,:} '/' denomIon{1,:}],'fontsize',12)
set(gca,'fontsize',12,'linewidth',1)

figure(8)
clf
mesh(x,y,transpose(matrixcompositions));
xlabel('N_m_i_n','fontsize',12)
ylabel('d_m_a_x','fontsize',12)
zlabel(['Matrix composition' numerIon{1,:} '/' denomIon{1,:}],'fontsize',12)
set(gca,'fontsize',12,'linewidth',1)
    
figure(9)
clf
mesh(x,y,transpose(compositionaca));
xlabel('N_m_i_n','fontsize',12)
ylabel('d_m_a_x','fontsize',12)
zlabel(['Cluster composition overall not averaged' numerIon{1,:} '/' denomIon{1,:}],'fontsize',12)
set(gca,'fontsize',12,'linewidth',1)

end

