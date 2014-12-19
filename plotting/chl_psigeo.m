%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     PLOTA VETORES DE VELOCIDADE GEOSTROFICA
%            ORIUNDAS DO CTD E ADCP
%   SOBRE UM MAPA DE CLOROFILA PARA A MESMA DATA
%
%        Rafael Soutelino - Mestrado IOUSP
%                  28/07/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% definindo os limites da area
lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');


% carregando arquivos com os dados de clorofila

load ../clorofila.mat
xchl = lon;
ychl = lat;
clear lon lat

[xchl,ychl] = meshgrid(xchl,ychl);
[xchl,ychl] = m_ll2xy(xchl,ychl,'patch');


% carregando arquivos com os dados de u e v geostroficos

load ../../hidrografia/mat/uv_psi_AO_OEII_1m.mat

load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

clear li lj xg yg seagrid_leste2 grid

% carregando arquivo com os dados de u e v oriundos do ADCP

load ../../adcp/mat/adcp_filt_AO.mat;

% montando mascara da isobata de 100 m

load ../../common/etopo2_leste.mat;
[c,h] = contour(xb,yb,zb,[-100 -100]);
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
xcont=xcont';
ycont=ycont';
close(1); 
xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

% plotando 

fc = 1;

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
stop
imagesc(xchl(1,:),ychl(:,1),chl,[-1.5 0])
axis('xy','equal','tight');hold on
m_quiver(xgi,ygi,ugo*fc,vgo*fc,0,'w');
fill(xf,yf,[0.8 0.8 0.8]);
m_usercoast('../../common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Clorofila + Velocidades Geostroficas - OEII','fontsize',10,'fontweight','bold');
drawnow
   print(1,'-depsc',['../figuras/clorofila_vgeo']);
   eval(['!epstopdf ../figuras/clorofila_vgeo.eps'])


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

imagesc(xchl(1,:),ychl(:,1),chl,[-1.5 0])
axis('xy','equal','tight');hold on
m_quiver(xgi,ygi,Uadcp*fc,Vadcp*fc,0,'w');
fill(xf,yf,[0.8 0.8 0.8]);
m_usercoast('../../common/costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('Clorofila + Velocidades ADCP - OEII','fontsize',10,'fontweight','bold');
drawnow
   print(2,'-depsc',['../figuras/clorofila_adcp']);
   eval(['!epstopdf ../figuras/clorofila_adcp.eps'])



















