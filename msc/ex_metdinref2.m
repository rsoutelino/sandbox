%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FUNCAO DE CORRENTE OBSERVADA PARA OS DADOS DO PRO-ABROLHOS I
%       A PARTIR DOS DADOS PROCESSADOS PELO CODAS - EXPERIMENTO 1
%            Rafael Soutelino - Mestrado IOUSP - ago/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% configuracoes ----------------------------------------------------
alis = 'n';
jan = 11;
cruz = 'oeii';
tr = 'oeii';
cc = 'n';
maperr = 'n';
%-------------------------------------------------------------------

% declarando projecao para trabalhar com m_map
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../proc/common/etopo2_leste.mat;

% preparando para a analise objetiva
load ../proc/common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

l1 = max(grid(:,1));
l2 = max(grid(:,2));

xgi = reshape(xg,l2,l1);
ygi = reshape(yg,l2,l1);

lc = 1; 
E = 0.02;

% criando u e v idealizados para dar exemplo do met din referenciado:
load coord.mat

dir1 = [220 220 220 220 220 180 170 90 180 250 180 100 180 240 180 160 220 220 220 220 220 180 170 90 180 250 180 100 180 240 180 160];
mod1 = [1.5 1.5 1.5 1.5 1.5 1.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.7 0.7 0.7 0.7 0.7 0.7 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];

dir2 = [40 40 40 40 40 0 170 90 180 250 180 100 180 240 180 160 40 40 40 40 40 0 170 90 180 250 180 100 180 240 180 160];
mod2 = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.2 0.2 0.2 0.2 0.2 0.2 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25];

[u1,v1] = intdir2uv(mod1,dir1,0,0);
[u2,v2] = intdir2uv(mod2,dir2,0,0);

uo1 = scaloa(xg,yg,lon,lat,u1,lc,E);
uo2 = scaloa(xg,yg,lon,lat,u2,lc,E);
vo1 = scaloa(xg,yg,lon,lat,v1,lc,E);
vo2 = scaloa(xg,yg,lon,lat,v2,lc,E);

uo1 = reshape(uo1,l2,l1);
uo2 = reshape(uo2,l2,l1);
vo1 = reshape(vo1,l2,l1);
vo2 = reshape(vo2,l2,l1);

uo3 = uo1 + uo2;
vo3 = vo1 + vo2;

s=0.3;

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');
hold on
m_quiver(xgi(1:1:end),ygi(1:1:end),uo1(1:1:end)*s,vo1(1:1:end)*s,0,'b');
[c,h]=m_contourf(xb,yb,zb,[-100 -100],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-100 -100],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('../proc/common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
title('v_{(150/10)}','fontsize',12,'fontweight','bold');
drawnow
print -depsc ../tex/figuras/ex_metdinref2_1.eps
!epstopdf ../tex/figuras/ex_metdinref2_1.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');
hold on
m_quiver(xgi(1:1:end),ygi(1:1:end),uo2(1:1:end)*s,vo2(1:1:end)*s,0,'b');
[c,h]=m_contourf(xb,yb,zb,[-100 -100],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-100 -100],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('../proc/common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
title('v_{obs(150)}','fontsize',12,'fontweight','bold');
drawnow
print -depsc ../tex/figuras/ex_metdinref2_2.eps
!epstopdf ../tex/figuras/ex_metdinref2_2.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');
hold on
m_quiver(xgi(1:1:end),ygi(1:1:end),uo3(1:1:end)*s,vo3(1:1:end)*s,0,'b');
[c,h]=m_contourf(xb,yb,zb,[-100 -100],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-100 -100],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('../proc/common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
title('v_{tot(10m)}','fontsize',12,'fontweight','bold');
drawnow
print -depsc ../tex/figuras/ex_metdinref2_3.eps
!epstopdf ../tex/figuras/ex_metdinref2_3.eps