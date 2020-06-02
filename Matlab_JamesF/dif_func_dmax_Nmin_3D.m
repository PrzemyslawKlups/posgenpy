%This scrip is to be use in conjuction with cluster sweep parameter
%selection. It highlights relevant regions of graph as described in the
%code.

%dinitial = 0.0;
%dstep =0.05;
%dfinal = 2;
%Nminmin = 2;
%Nminmax =20;
Nminnp = Nminmax - Nminmin +1;
dmaxnp = round((dfinal-dinitial+dstep)/dstep);%number of points for dmax%    
dif = zeros(Nminnp,dmaxnp);
numreal = zeros(Nminnp,dmaxnp);
for i=1:1:Nminnp
    name =['Nmin' num2str(Nminmin -1 + i)]; %vaule of dmax adding zeros in for%
    b = Ncdifhighp_Struct.(name);%extract relevant array from structure%
    a = b(:,1); %extract frequency column for that array%
    dif(i,1:size(b(:,1)))= a; %difference into matrix of zeros, maxb bit makes dimensions match up%
end
[x,y]= meshgrid(Nminmin:Nminmax,dinitial:dstep:dfinal); %mesh of possible dmax and Nmin values values%
figure(10)
mesh(x,y,dif');
xlabel('Nmin','fontsize',12)
ylabel('dmax','fontsize',12)
zlabel('difference in number clusters','fontsize',12)
set(gca,'fontsize',12,'linewidth',1)

for i=1:1:Nminnp
    name =['Nmin' num2str(Nminmin -1 + i)]; %vaule of dmax adding zeros in for%
    b = Nchighp_Struct.(name);%extract relevant array from structure%
    a = b(:,1); %extract frequency column for that array%
    numreal(i,1:size(b(:,1)))= a; %difference into matrix of zeros, maxb bit makes dimensions match up%
end
[x,y]= meshgrid(Nminmin:Nminmax,dinitial:dstep:dfinal); %mesh of possible dmax and Nmin values values%
figure(11)
mesh(x,y,numreal');
xlabel('Nmin','fontsize',12)
ylabel('dmax','fontsize',12)
zlabel('number clusters','fontsize',12)
set(gca,'fontsize',12,'linewidth',1)
ylim([0,2])

figure(3)
hold on
surf(x,y,numreal','FaceColor','r'); %highlights the valid section of the graph%
hold off


