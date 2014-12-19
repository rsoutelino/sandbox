%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              PROBLEMA  3                %%%
%%%               WESTRAX 2                 %%%
%%%     Rafael Soutelino - julho / 2006     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc

%%% LOOPING PARA LEITURA DOS DADOS

lon=[]; lat=[]; nest=[]; cont=1; Tsec=[]; Ssec=[];

for i=[35 36 37 39 40 41 42] 

    eval(['load ../ctd/w2cpz0',num2str(i),'.pro;']);
    eval(['data = w2cpz0',num2str(i),';'])
    
    p=data(:,1);
    f=find(p==1500);
    if isempty(f)==0;
       lon=[lon data(1,6)*-1];
       lat=[lat data(1,5)];
       eval(['p',num2str(i),' = data(5:3002,1);'])
       eval(['u',num2str(i),' = data(5:3002,8);'])
       eval(['v',num2str(i),' = data(5:3002,9);'])
       nest=[nest i];
    end

    eval(['clear w2cpz0',num2str(i)])    
end

clear data   

Usec=[u42 u41 u40 u39 u37 u36 u35];
Vsec=[v42 v41 v40 v39 v37 v36 v35];
psec=-p35;

dist=sw_dist(lat,lon,'km');
dist=[0 cumsum(dist)];


%%% decimando N2 e p

U=[]; V=[]; p=[];

for i=10:10:length(psec)
    U(i/10,:) = Usec(i,:);
    V(i/10,:) = Vsec(i,:);
    p(i/10,:) = psec(i,:);
end

lU=[-90:10:100];


figure

[c,h]=contourf(dist,p,U,lU); hold on; colorbar
%  shading flat;
axis([0 max(dist) -1500 -10])
%  clabel(c,h)
title('Secao Vertical de Velocidades - Pegasus','fontweight','bold')
xlabel('Distancia da Costa')
ylabel('Profundidade (m)')
%  set(gca,'plotboxaspectratio',[.5 1 .5])
set(gcf,'color','w')

print -depsc pegasus.eps
!epstopdf pegasus.eps

clc







