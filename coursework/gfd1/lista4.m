clear all;close all;clc

t = -20:0.01:20;

f0 = 1;
r = 0.1;
x0 = 1;
y0 = 1;

x = (x0./r).*f0.*exp(r.*t).*(sin(f0.*t) + (r/f0).*cos(f0.*t));
y = y0.*exp(r.*t).*((r/f0).*sin(f0.*t) - cos(f0.*t));

figure(1)

set(gcf,'color','w')

subplot(311)
plot(fliplr(t),x,'linewidth',2)
axis([-20 20 -60 60])
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])
ylabel('x','fontweight','bold')
title('Componente zonal da trajetoria','fontweight','bold')
pbaspect([1 .5 1])

subplot(312)
plot(fliplr(t),y,'r','linewidth',2)
axis([-20 20 -7 6])
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])
ylabel('y','fontweight','bold')
xlabel('t','fontweight','bold')
title('Componente meridional da trajetoria','fontweight','bold')
pbaspect([1 .5 1])

subplot(313)
plot(x,y,'k','linewidth',1.5)
axis([-70 70 -8 8])
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])
ylabel('y','fontweight','bold')
xlabel('x','fontweight','bold')
title('Trajetoria no plano horizontal','fontweight','bold')
pbaspect([1 .5 1])

print -depsc lista4_exe6.eps
!epstopdf lista4_exe6.eps
!rm -rf lista4_exe6.eps




