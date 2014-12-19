%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%    Lista Pratica # 2 - Problemas 1 e 2                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
format long g 

lons=[]; lats=[]; 

xx=[35:37 39:42];lat=5;
lx=length(xx);

P=1:4000;
TT=NaN*ones(length(P),lx);
SS=NaN*ones(length(P),lx);

for i=1:lx
    
    eval(['load ../ctd/w2cpz0',num2str(xx(i)),'.pro']);
    eval(['est=w2cpz0',num2str(xx(i)),';']);
    eval(['clear w2cpz0',num2str(xx(i))]);
    lons=[lons est(1,6)*-1];
    lats=[lats est(1,5)];

p=est(2:end,1);
f=find(p>=P(1) & p<=P(end));

TT(1:length(f),i)=est(f,3);
SS(1:length(f),i)=est(f,5);
u(1:length(f),i)=est(f,8);
v(1:length(f),i)=est(f,9);

end

u=fliplr(u);
v=fliplr(v);

dist=sw_dist(lats,lons,'km');
dist=[0 cumsum(dist)];

Tm=nanmean(TT')';
%eval(['Tms=weim(41,''hann'',Tm);']);
Sm=nanmean(SS')';
%eval(['Sms=weim(41,''hann'',Sm);']);

[bfrq,vort,p_ave] = sw_bfrq(Sm,Tm,P',[5]);
bfrq=abs(bfrq);
eval(['bfrqs=weim(21,''hann'',bfrq);']);

pbin=[10:10:max(p_ave)+.5]';
ln=length(pbin);

     for kk=1:ln
         f=find(p_ave>=pbin(kk)-4.5 & p_ave<=pbin(kk)+4.5);
         bfm(kk)=mean(bfrqs(f));
     end

n2=bfm';


%***************************************************************
%           THE EIGENFUNCTION CALCULATION
%***************************************************************
% COMPUTING THE 10-METER AVG N2 PROFILE:

global f;
global xi0;
global po;
global n2f2; % Usado p/ calcular o estiramento vertical


po=pbin(1:end-1)+5;

f0=14.585e-5*sin(lat*pi/180);
beta=14.585e-5*cos(lat*pi/180)/6400000;
f2=f0*f0;


% EIGENMODE COMPUTATION

n2f2=n2/f2;
rd=[197.7 126.2 74.1 54.6 44.4];
eigm=zeros(length(po),5);

for i=1:5,
  r(i)=fzero('vmode2',rd(i));
  eigm(:,i)=f';
end
r
stop
r0=1e-3*sqrt(9.81*4000)/f0

h0=((r0*1e3)^2*f2)/9.81;
h1=((r(1)*1e3)^2*f2)/9.81;
h2=((r(2)*1e3)^2*f2)/9.81;
h3=((r(3)*1e3)^2*f2)/9.81;
h4=((r(4)*1e3)^2*f2)/9.81;
h5=((r(5)*1e3)^2*f2)/9.81;


r=r';
F1=eigm(:,1);
F2=eigm(:,2);
F3=eigm(:,3);
F4=eigm(:,4);
F5=eigm(:,5);
F0=ones(length(F1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PLOTTING                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CONVERSAO DE N2 DE RAD2/SEC2 TO CPH %%%

bv1=sqrt(n2)*3600/(2*pi);
z1=-pbin;
zer = [0 -10;0 -4000];

%subplot(211)
figure
plot(zer(:,1),zer(:,2), 'k-')
axis([0 15 -4000 0])
axis('square')
hold on
plot(bv1,z1,'k','linewidth',2)
NAME=['WX2 Radial 44^oW - Frequencia de Estratificacao N(z)'];
title(NAME,'fontsize',12,'fontweight','bold')
xlabel('N [cph]','FontSize',10)
ylabel('Profundidade [m]','FontSize',10)
set(gca,'YTick',[-4000:500:0],'YTickLabel',[4000:-500:0]);
set(gcf,'color','w')
hold off



figure
subplot(131)
plot(Tm,-P,'r','linewidth',2)
xlabel('Temperatura [^oC]','FontSize',10)
ylabel('Profundidade [m]','FontSize',10)
NAME=['Temperatura'];
title(NAME,'fontsize',12,'fontweight','bold')
axis([0 28 -4000 0])
set(gcf,'color','w')
set(gca,'plotboxaspectratio',[.5 1 .5])

subplot(132)
plot(Sm,-P,'g','linewidth',2)
xlabel('Salinidade','FontSize',10)
%  ylabel('Profundidade [m]','FontSize',10)
NAME=['Salinidade'];
title(NAME,'fontsize',12,'fontweight','bold')
axis([34 36.8 -4000 0])
set(gca,'plotboxaspectratio',[.5 1 .5])


subplot(133)
plot(bv1,z1,'b','linewidth',2);hold on
NAME=['Frequencia de Estratificacao N(z)'];
title(NAME,'fontsize',12,'fontweight','bold')
xlabel('N [cph]','FontSize',10)
%  ylabel('Profundidade [m]','FontSize',10)
set(gca,'YTick',[-4000:500:0],'YTickLabel',[4000:-500:0]);
axis([0 15 -4000 0])
set(gca,'plotboxaspectratio',[.5 1 .5])

print -depsc n2.eps
!epstopdf n2.eps

% PLOT THE PRESSURE EIGENMODES

%subplot(212)
figure
set(gcf,'color','w')
plot(zer(:,1),zer(:,2),'-','Color',[0.753 0.753 0.753])
axis([-5 6 -4000 0])
axis('square')
hold on
plot(F0,-po,'r-.','linewidth',2)
plot([3.5 4],[-2700 -2700],'r-.','linewidth',2);text(4.25,-2700,'Zero');
plot(F1,-po,'b-','linewidth',2)
plot([3.5 4],[-2900 -2900],'b-','linewidth',2);text(4.25,-2900,'Primeiro');
plot(F2,-po,'g--','linewidth',2)
plot([3.5 4],[-3100 -3100],'g--','linewidth',2);text(4.25,-3100,'Segundo');
plot(F3,-po,'m-')
plot([3.5 4],[-3300 -3300],'m-');text(4.25,-3300,'Terceiro');
plot(F4,-po,'y--')
plot([3.5 4],[-3500 -3500],'y--');text(4.25,-3500,'Quarto');
plot(F5,-po,'c-.')
plot([3.5 4],[-3700 -3700],'c-.');text(4.25,-3700,'Quinto');
title('Estrutura dos Seis Primeiros Modos Dinamicos Verticais','fontsize',12,'fontweight','bold') 
xlabel('Amplitude dos Modos','fontsize',10)
ylabel('Profundidade [m]','fontsize',10)
set(gca,'YTick',[-4000:500:0],'YTickLabel',[4000:-500:0]);
hold off

print -depsc modos.eps
!epstopdf modos.eps

%%% SALVANDO OS RAIOS DE DEFORMACAO EM ARQUIVO TEXTO %%%

RR0=['Rd_0='];RR1=['Rd_1='];RR2=['Rd_2='];RR3=['Rd_3='];RR4=['Rd_4='];RR5=['Rd_5='];
RRt=[RR0;RR1;RR2;RR3;RR4;RR5];
RRs=[r0;r(1);r(2);r(3);r(4);r(5)];
RRs=num2str(RRs);
pl=['\n';'\n';'\n';'\n';'\n';'\n'];
RR=[RRt RRs pl];
fid=fopen('raios.txt','w');
fprintf(fid,RR');
fclose(fid);

%%% SALVANDO AS PROFUNDIDADES EQUIVALENTES EM ARQUIVO TEXTO %%%

hn=[h0;h1;h2;h3;h4;h5];
hr=num2str(hn);
hR=['h0=';'h1=';'h2=';'h3=';'h4=';'h5='];
HR=[hR hr pl];
fid=fopen('profequi.txt','w');
fprintf(fid,HR');
fclose(fid);


save F F0 F1 F2 F3 F4 F5 po

%%% SECOES

% decimando 

U=[]; V=[]; p=[];

for i=10:10:(length(P)-1)
    U(i/10,:) = u(i,:);
    V(i/10,:) = v(i,:);
    p(i/10,:) = P(i);
end


lu=[-90:10:100];

figure

[c,h]=contourf(dist,-p,U,lu); hold on; colorbar('horiz')
shading flat;
axis([0 max(dist) -1500 -10])
clabel(c,h)
title('Secao Vertical de Velocidades - Pegasus','fontweight','bold')
xlabel('Distancia da Costa')
ylabel('Profundidade (m)')
%  set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

print -depsc pegasus.eps
!epstopdf pegasus.eps

%%% calculando as amplitudes Ui

U=U/100; % passando para m/s

for i=0:5
   for k=1:length(xx)
       eval(['Amod',num2str(i),'(1,k) = 10/4000*sum(U(:,k).*(F',num2str(i),'));'])
   end
end

for i=0:5
   for k=1:length(xx)
       eval(['Umod',num2str(i),'(:,k) = Amod',num2str(i),'(1,k).*F',num2str(i),';'])
   end
end


figure

[c,h]=contourf(dist,-p,Umod0*100+Umod1*100,lu); hold on; colorbar('horiz')
shading flat;
axis([0 max(dist) -1500 -10])
clabel(c,h)
title('2 primeiros modos','fontweight','bold')
xlabel('Distancia da Costa')
ylabel('Profundidade (m)')
%  set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

print -depsc 2modos.eps
!epstopdf 2modos.eps

figure

[c,h]=contourf(dist,-p,Umod0*100+Umod1*100+Umod2*100,lu); hold on; colorbar('horiz')
shading flat;
axis([0 max(dist) -1500 -10])
clabel(c,h)
title('3 primeiros modos','fontweight','bold')
xlabel('Distancia da Costa')
ylabel('Profundidade (m)')
%  set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

print -depsc 3modos.eps
!epstopdf 3modos.eps


figure

[c,h]=contourf(dist,-p,Umod0*100+Umod1*100+Umod2*100+Umod3*100+Umod4*100+Umod5*100,lu); hold on; colorbar('horiz')
shading flat;
axis([0 max(dist) -1500 -10])
clabel(c,h)
title('6 primeiros modos','fontweight','bold')
xlabel('Distancia da Costa')
ylabel('Profundidade (m)')
%  set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

print -depsc 6modos.eps
!epstopdf 6modos.eps


clc
























