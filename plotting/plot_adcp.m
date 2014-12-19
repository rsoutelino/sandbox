% plotagem rapida ADCP oceano leste 2
clear all;close all;clc

u = []; v = []; lon = []; lat = [];


lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');

for k=1:54
    eval(['load adcp',num2str(k),'.mat;'])
    lon = [lon AnFLonDeg'];
    lat = [lat AnFLatDeg'];
    u = [u SerEmmpersec'];
    v = [v SerNmmpersec'];
end

f = find(u == -32768);
u(f) = nan;
v(f) = nan;

u = u(1,1:10:end);
v = v(1,1:10:end);
lon = lon(1:10:end);
lat = lat(1:10:end);

f = find(isnan(u)==0);
u = u(f); v = v(f); lon = lon(f); lat = lat(f);
u = weim(31,'hann',u);
v = weim(31,'hann',v);

load ../../ctd/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

xg=xg';yg=yg';
corrlen = 4; 
err = 0.005; 

load ../../ctd/etopo2_leste2.mat;
n=1
% buscando valor de lon e lat das isobatas de ?? m
if n <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-n -n]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
xcont = weim(31,'hann',xcont);
ycont = weim(31,'hann',ycont);
%  xcont=xcont';
%  ycont=ycont';

close(1); 

xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

% enxertando valores da isobata nos vetores
lon = [lon xcont]; lat = [lat ycont];
u = [u zeros(size(xcont))];
v = [v zeros(size(xcont))];

[psiob]=vectoa(xg,yg,lon,lat,u,v,corrlen,err,0);

psiob=reshape(psiob,li,lj);
xg=reshape(xg,li,lj);
yg=reshape(yg,li,lj);

[uo,vo]=psi2uv(xgi,ygi,psiob);

fc = 405; % fator de escala dos vetores

figure(4)
set(gcf,'color','w')
m_contourf(xg,yg,psiob,[min(min(psiob)):10:max(max(psiob))]);hold on;shading flat;
%  m_contour(xgi,ygi,psigo,[min(min(psigo)):1000:max(max(psigo))],'b')
m_quiver(xg,yg,uo*fc,vo*fc,0,'k')
fill(xf,yf,[0.8 0.8 0.8]);
%  plot(xf,yf,'color',[.8 .8 .8])
%  m_plot(lonctd,latctd,'w.','markersize',6)
m_usercoast('costa_leste.mat','patch',[0 0 0])
vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
m_quiver(-40,-11,vmax*fc/100,0,0,'w')
text = [num2str(vmax),' cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',8,'fontweight','bold')
%  text = [nn,' m'];
%  m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
title('\Psi Observado [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1)/1.2 pos(2) pos(3)/2.5 pos(4)])
   print -depsc teste.eps
drawnow
