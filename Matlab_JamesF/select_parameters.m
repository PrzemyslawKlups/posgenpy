function [Ncdifhighp_Struct,Nchighp_Struct,percent_real_Struct] = select_parameters(dinitial,dstep,dfinal,Nminmin,Nminmax,sizedist_Struct,sizedistRL_Struct,p)
%select parameters takes the size distribution files from a cluster sweep
%and processes
%Inputs
%dinitial, dstep, dfinal, - values of dmax to plot (NB can be a subset of all the values used in cluster sweep
%Nminmin, Nminmax - max and min values of Nmin to plot, step size is 1
%sizedist_Struct and sizedist_Struct_RL - dataStructures containing all the information from saved files, made by form_sizedist_Struct & form_sizedist_Struct_RL functions
% p - threshold value, fraction of real clusters above which dmax and Nmin are excpted
%

%Outputs
% Ncdifhighp_Struct - Diference between real and random, but only in the region where precent real > p*100
% Nchighp_Struct - Number of clusters, but only in the region where precent real > p*100
% percent_real_Struct - datastructure containing percentage of real clusters, (Nreal - Nrand)/Nreal*100

%NB p is a fraction, percent_real_Struct which it is used to threshold is a
%percentage - sorry

%Plots
%figure(4) - percent_real_Struct, done as multiple 2D lines not 3D
%figure (5)and figure(6) - plot differnece in number of clusters between real and random, 
% somewhat less useful. I did also consider using no of atoms clustered
% instead of number of clusters and this relates to the work comparing
% those. Can comment out if not required

dlargest = ['dmax' num2str(dfinal*100)];
maxcluster = size(sizedist_Struct.(dlargest),1);%largest cluster (for all values of dmax)%
dmaxnp = round((dfinal-dinitial+dstep)/dstep);%number of points for dmax%    
f= zeros(dmaxnp,maxcluster);
maxclusterRL = size(sizedistRL_Struct.(dlargest),1);%largest cluster (for all values of dmax)%
dmaxnpRL = round((dfinal-dinitial+dstep)/dstep);%number of points for dmax%  
fRL = zeros(dmaxnpRL,maxclusterRL);
for i=1:1:dmaxnp
    name =['dmax' num2str(dinitial*100+i*dstep*100-(dstep*100))]; %vaule of dmax adding zeros in for%
    b = sizedist_Struct.(name);%extract relevant array from structure%
    a =b(:,2); %extract frequency column for that array%
    f(i,1:max(b(:,1)))= a; %frequencies into matrix of zeros, maxb bit makes dimensions match up%
end
for i=1:1:dmaxnpRL
    name =['dmax' num2str(dinitial*100+i*dstep*100-(dstep*100))]; %vaule of dmax adding zeros in for%
    b = sizedistRL_Struct.(name);%extract relevant array from structure%
    a =b(:,2); %extract frequency column for that array%
    fRL(i,1:max(b(:,1)))= a; %frequencies into matrix of zeros, maxb bit makes dimensions match up%
end
for i=Nminmin:Nminmax %range of Nmin values plotting%
    roif = f(:,i:size(f,2)); %region of interest in terms on clipping out small N%
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    Nreal = sum(roif,2); %total no. clusters for each dmax value%
    roifRL = fRL(:,i:size(fRL,2));
    Nrand = sum(roifRL,2);
    percent_real = ((Nreal-Nrand)./Nreal).*100; %selection criteria%
    Nclusterdif = Nreal-Nrand; %difference in Number of clusters%
    %form data structure%
    percent_real_Struct.(fieldNmin)= percent_real;
    Nclusterdif_Struct.(fieldNmin) = Nclusterdif;
    Ncluster_real_Struct.(fieldNmin) = Nreal;
end
x = dinitial:dstep:dfinal;
figure(4)
clf
for i=Nminmin:Nminmax
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    hold on
    figure(4)
    plot(x,percent_real_Struct.(fieldNmin),'Color',[sin(pi/2*((i-Nminmin)/(Nminmax-Nminmin))),sin(pi*((i-Nminmin)/(Nminmax-Nminmin))),cos(pi/2*((i-Nminmin)/(Nminmax-Nminmin)))]);
    hold off
end
figure(4)
xlabel('dmax','fontsize',12)
ylabel('Percentage of clusters not expected','fontsize',12)
leg = {1, Nminmax - Nminmin + 1};
for i = 1:Nminmax - Nminmin + 1
    leg{1,i} = ['Nmin = ' num2str(i + Nminmin -1)];
end
legend(leg);

%now start filtering based on figure 4%
clear Ncdifhighp_Struct
for i=Nminmin:Nminmax
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    a = percent_real_Struct.(fieldNmin);
    b = Nclusterdif_Struct.(fieldNmin);
    for j=1:size(a,1)
        if a(j,1) < p*100
            b(j,1)=NaN;
        end
        Ncdifhighp_Struct.(fieldNmin)= b;%region of difference in no. clusters where cluster percent is above threshold%       
    end
    
end
x = dinitial:dstep:dfinal;
figure(5)
clf
for i=Nminmin:Nminmax
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    figure(5)
    hold on
    plot(x,Nclusterdif_Struct.(fieldNmin),'--','Color',[sin(pi/2*((i-Nminmin)/(Nminmax-Nminmin))),sin(pi*((i-Nminmin)/(Nminmax-Nminmin))),cos(pi/2*((i-Nminmin)/(Nminmax-Nminmin)))]);
    plot(x,Ncdifhighp_Struct.(fieldNmin),'Color',[sin(pi/2*((i-Nminmin)/(Nminmax-Nminmin))),sin(pi*((i-Nminmin)/(Nminmax-Nminmin))),cos(pi/2*((i-Nminmin)/(Nminmax-Nminmin)))]);
    hold off
end
figure(5)
xlabel('dmax','fontsize',12)
ylabel('Difference in no. clusters and expected','fontsize',12)
figure(6)
clf
for i=Nminmin:Nminmax
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    figure(6)
    hold on
    plot(x,Ncdifhighp_Struct.(fieldNmin),'Color',[sin(pi/2*((i-Nminmin)/(Nminmax-Nminmin))),sin(pi*((i-Nminmin)/(Nminmax-Nminmin))),cos(pi/2*((i-Nminmin)/(Nminmax-Nminmin)))]);
    hold off
end
figure(6)
xlabel('dmax','fontsize',12)
ylabel('Difference in no. clusters and expected','fontsize',12)
legend(leg,'Location','northwest');
%Number of clusters, when I've filtered out unacceptable values of dmax and
%Nmin%
clear Nchighp_Struct
for i=Nminmin:Nminmax
    fieldNmin = ['Nmin' num2str(i)];%field name for structure%
    a = percent_real_Struct.(fieldNmin);
    b = Ncluster_real_Struct.(fieldNmin);
    for j=1:size(a,1)
        if a(j,1) < p*100
            b(j,1)=NaN;
        end
        Nchighp_Struct.(fieldNmin)= b;%region of difference in no. clusters where cluster percent is above threshold%       
    end
    
end
end

