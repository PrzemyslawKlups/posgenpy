function [sizedist_Struct] = form_sizedist_Struct(dinitial,dstep,dfinal,inputfolder,plot_f)
%form_sizedist_Struct reads the size dist files produced by a cluster sweep
%and form a data structure with the data.
%fieldnames in the structure are of the form dmax0.55
%each field is something by 2 double with cluster size and number of that
%size. Includes 0,1 in each structure
%
%inputfolder can just be output folder from ClusterSweep%
clear sizedist_Struct %clear structures as data doesn't overwrite%
Nmin = 2; %writting these here now even though intention is for this to not be altered at this stage as it makes editing easier%
for i = dinitial:dstep:dfinal
    sizedistfile = [inputfolder '\sizedist-dmax-' num2str(i) '.txt'];
    name =['dmax' num2str(i*100)];%name of field in dataStruct, name is 100x dmax as field name can't include .%
    fileID = fopen(sizedistfile);
    C = textscan(fileID,'%s','Delimiter', '');
    fclose(fileID);
    sz = size(C{1,1},1);
    if sz>1
    sizedist_Struct.(name) = dlmread(sizedistfile,'',1,0);%combining all the data into datastruct%
    else
    sizedist_Struct.(name) = [2,0];%this is because it doesn't like reading files where no clusters were found%
    end
end 
for j=dinitial:dstep:dfinal %add zeros in where there are no cluster%
    name =['dmax' num2str(j*100)]; %vaule of dmax adding zeros in for%
    a = sizedist_Struct.(name); %saving typing dataStruct everytime%
    a0 = zeros(max(a(:,1)),size(a,2));%matrix of zeros of correct dimensions%
    a0(a(:,1),:)= a;%adding in rows with zeros%
    a0(:,1)=Nmin-1:1:max(a(:,1));%making clusters size values non-zero%
    sizedist_Struct.(name)=a0;%reforming data Structure%
end

if strcmp(plot_f,'Yes') == 1
    dlargest = ['dmax' num2str(dfinal*100)];
    maxcluster = size(sizedist_Struct.(dlargest),1);%largest cluster (for all values of dmax)%
    dmaxnp = round((dfinal-dinitial+dstep)/dstep);%number of points for dmax%
    dmaxrng = dinitial:dstep:dfinal;% dmax range for plot%
    clusterszrng = Nmin-1:1:maxcluster;%linearly spaces cluster size values%
    [x,y]= meshgrid(clusterszrng,dmaxrng); %mesh of possible dmax and cluster size values%
    f= zeros(dmaxnp,maxcluster);
    for i=1:1:dmaxnp
        name =['dmax' num2str(dinitial*100+i*dstep*100-(dstep*100))]; %vaule of dmax adding zeros in for%
        b = sizedist_Struct.(name);%extract relevant array from structure%
        c =b(:,2); %extract frequency column for that array%
        f(i,1:max(b(:,1)))= c; %frequencies into matrix of zeros, maxb bit makes dimensions match up%
    end 
    figure(1)
    clf
    mesh(x,y,f);
    xlim([0,maxcluster])
    ylim([dinitial,dfinal])
    xlabel('cluster size','fontsize',12)
    ylabel('dmax','fontsize',12)
    zlabel('frequency real','fontsize',12)
    set(gca,'fontsize',12,'linewidth',1)
end
end

