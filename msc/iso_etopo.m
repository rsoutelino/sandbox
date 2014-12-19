%%% extrai coordenadas de isobatas a partir
%%% do banco ETOPO2

clear all;close all;clc

load etopo2_E_751x901.xyz;

x=etopo2_E_751x901(:,1);
y=etopo2_E_751x901(:,2);
z=etopo2_E_751x901(:,3);

x = reshape(x,751,901);
y = reshape(y,751,901);
z = reshape(z,751,901);

contour(x,y,z); close(1);
%  z=smoo2(z,-9999,9,0.1);
[c,h]=contour(x,y,z,[-1000 -1000]);
x1000=get(h(1),'xdata'); y1000=get(h(1),'ydata');

save iso1000.mat x1000 y1000
close(1)

% limitando uma area qualquer

xmin = -42;
xmax = -35;
ymin = -23;
ymax = -14;

fx = find(x(:,1) >= xmin & x(:,1) <= xmax);
fy = find(y(1,:) >= ymin & y(1,:) <= ymax);

x = x(fx,fy);
y = y(fx,fy);
z = z(fx,fy);

f = find(z > 0);
z(f) = 0; 

load ../proc/common/costa_leste.mat

xlc = ncst(:,1);
ylc = ncst(:,2);

f2 = find(ylc >= ymin-1 & ylc <= ymax);

xlc = xlc(f2);
ylc = ylc(f2);

figure(1)
set(gcf,'color','w')
surf(x,y,z); shading interp;hold on
plot(xlc,ylc,'k.','markersize',0.2)
xlabel('Longitude')
ylabel('Latitude')
zlabel('Profundidade')
cc=colorbar('horiz');
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2) pos(3) pos(4)/2])
print -dpng batimetria_abrolhos.png 

