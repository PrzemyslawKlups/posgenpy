function [varargout] = Total_clusters_func_dmax_Nmin(dinitial, dstep, dfinal, Nminmin, Nminmax, sizedist_Struct)
%Total_cluster_func_dmax_Nmin Summary - takes a size distribution structure
%and creates a 3D plot of total number of clusters as a function od dmax
%and Nmin
% Usual dinitial, dstep, dfinal inputs as well as a user choice on the size
% of clusters they what to plot - NB clusters larger than Nminmax are still
% counted towards the total.

%NB Z axis label is hardcoded for Total number of clusters real, edit as
%appropriate%
%NB also harded coded to plot figure 3, also edit as appropriate%

NumberCm = zeros(round((dfinal- dinitial + dstep)/dstep),Nminmax-Nminmin+1);%Number Clusters matrix%
%matrix of zeros to be filled by a value of total no. cluster so the given dmax and Nmin value%
dlargest = ['dmax' num2str(dfinal*100)];
maxcluster = size(sizedist_Struct.(dlargest),1);%largest cluster (for all values of dmax)%
dmaxnp = round((dfinal-dinitial+dstep)/dstep);%number of points for dmax%    
f= zeros(dmaxnp,maxcluster);
for i=1:1:dmaxnp
    name =['dmax' num2str(dinitial*100+i*dstep*100-(dstep*100))]; %vaule of dmax adding zeros in for%
    b = sizedist_Struct.(name);%extract relevant array from structure%
    c =b(:,2); %extract frequency column for that array%
    f(i,1:max(b(:,1)))= c; %frequencies into matrix of zeros, maxb bit makes dimensions match up%
end 
for i=Nminmin:Nminmax %range of Nmin values plotting%
    roif = f(:,i:size(f,2)); %region of interest in terms on clipping out small N%
    NumberC = sum(roif,2); %total no. clusters for each dmax value%
    %form matrix column all value of Nreal and Nrand for one value of Nmin structure%
    NumberCm(:,i - Nminmin + 1) = NumberC;
    
end
[x,y]= meshgrid(Nminmin:Nminmax,dinitial:dstep:dfinal); %mesh of possible dmax and Nmin values values%
figure(3)
clf
mesh(x,y,NumberCm);
xlabel('N_m_i_n','fontsize',12)
ylabel('d_m_a_x','fontsize',12)
zlabel('Total Number of clusters real','fontsize',12)
set(gca,'fontsize',12,'linewidth',1)

if nargout == 1
    varargout{1} = NumberCm; %Number Clusters matrix%
end
end

