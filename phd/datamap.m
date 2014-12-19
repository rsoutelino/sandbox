clear all; close all; clc;

lonlim = [-45 -30]; latlim = [-24 -9];
m_proj('mercator','lat',latlim,'lon',lonlim);
[zb,xb,yb] = m_tbase([lonlim latlim]);
m_gshhs_l('save','costa.mat');

m_usercoast('costa.mat','patch',[.9 .9 .9]); hold on
m_contour(xb,yb,zb,[-200 -1000],'k');
m_grid

load leste1/ctd/posicoes.dat; p = posicoes;
lon = p(:,3); lat = p(:,2);
m_plot(lon,lat,'.r')

load abrolhos1/ctd/posicoes.dat; p = posicoes;
lon = -[p(:,3) + p(:,4)/60]; lat = -[p(:,1) + p(:,2)/60];
m_plot(lon,lat,'.b')

load leste2/ctd/posicoes.dat; p = posicoes;
lon = -[p(:,4) + p(:,5)/60]; lat = -[p(:,2) + p(:,3)/60];
m_plot(lon,lat,'or')

load abrolhos2/ctd/posicoes.dat; p = posicoes;
lon = -[p(:,3) + p(:,4)/60]; lat = -[p(:,1) + p(:,2)/60];
m_plot(lon,lat,'ob')

load proabrolhos/ctd/posicoes.dat; p = posicoes;
lon = p(:,2); lat = p(:,3);
m_plot(lon,lat,'*k')

print -depsc datamap.eps

