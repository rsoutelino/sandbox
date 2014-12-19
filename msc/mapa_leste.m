clear all; close all; clc

load ../proc/common/etopo2_leste.mat;
load ../proc/common/costa_seagrid.mat;

f = find(lat >= -23 & lat <= -10);
lon = lon(f);
lat = lat(f);

zb = zb-40;
f2 = find(zb>0); zb(f2)=0;

load ../proc/common/floor_colormap2

figure(1)
set(gcf,'color','w')
surf(xb,yb,zb); shading interp; hold on
caxis([-9000 100])	
colormap(floor_colormap2)
p=plot(lon,lat,'k.','markersize',1);
%  set(p,'color',[.7 .7 .7])
%  p=plot(xb(f2),yb(f2),'k.');
%  set(p,'color',[.7 .7 .7])
% axis([-43 -30 -24 -10 -6000 1000])
axis([-43 -30 -24 -18.5 -6000 1000])
view([1 -2.5 15])
view([17.5 78])
h = light;
set(h,'Position',[10 -60 40000]);
material dull
text(-36,-17.7,'Banco Hot Spur','color','w','fontsize',10,'fontweight','bold')
text(-37.7,-16,'Banco Royal Charlote','color','w','fontsize',10,'fontweight','bold')
text(-39.5,-18.6,'Banco de Abrolhos','color','w','fontsize',10,'fontweight','bold')
text(-37.5,-21.5,'Cadeia Vitoria-Trindade','color','w','fontsize',10,'fontweight','bold')
text(-40,-13,100,'Salvador','color','w','fontsize',10,'fontweight','bold')
text(-42.5,-22.2,100,'C. Sao Tome','color','w','fontsize',10,'fontweight','bold')
text(-41,-20,100,'Vitoria','color','w','fontsize',10,'fontweight','bold')
xlabel('longitude [^\circ W]','fontweight','bold')
ylabel('latitude [^\circ S]','fontweight','bold')
zlabel('profundidade [m]','fontweight','bold')

print -depsc ../figuras/mapa_leste.eps
print -dpng ../figuras/mapa_leste.png
