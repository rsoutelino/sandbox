clear all;close all;clc

% criando perfil idealizado de velocidade

z1 = -2.3:0.001:1.5;
r = randn(size(z1));
r = r*0.1;
r = weim(91,'hann',r);

v = sin(z1);
v = v*0.5+0.2;
%  v = v+r;
v = fliplr(-v);

z = 1:length(z1);
z = z/max(z);
z = z*1000;

f = find(round(z)==150);
v2 = v - v(f(1)).*ones(size(v));

f2 = find(round(z)==1000);
v3 = v - v(f2(1)).*ones(size(v));


figure
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 10 4],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 10 4],...
        'ShareColors','off',...
        'Clipping','on');

subplot(131)
p = plot([0 0],[-10000 100],'k--');hold on
set(p,'color',[.7 .7 .7])
plot([-10 10],[-1000 -1000],'g')
plot(v3,-z,'linewidth',2);
axis([-1 1 -1000 0])
title('V_g --> p_0 = 1000 dbar','fontweight','bold')
ylabel('Profundidade [m]','fontweight','bold')
t = ['V_g(1000m) = 0'];
text(0.15,-100,t,'fontweight','bold','fontsize',8)

subplot(132)
p = plot([0 0],[-10000 100],'k--');hold on
set(p,'color',[.7 .7 .7])
plot([-10 10],[-150 -150],'g')
plot(v2,-z,'linewidth',2);
axis([-1 1 -1000 0])
title('V_g --> p_0 = 150 dbar','fontweight','bold')
xlabel('Velocidade [m s^{-1}]','fontweight','bold')
t = ['V_g(150m) = 0'];
text(0.15,-100,t,'fontweight','bold','fontsize',8)

subplot(133)
p = plot([0 0],[-10000 100],'k--');hold on
set(p,'color',[.7 .7 .7])
plot([-10 10],[-150 -150],'g')
plot(v,-z,'r','linewidth',2);
arrow([0 -150],[-0.6 -150],'length',1,'tipangle',7,'linewidth',1,'facecolor','k','edgecolor','k')
axis([-1 1 -1000 0])
title('V_g --> ref. ADCP','fontweight','bold')
t = ['V_{adcp}(150m) = -0.6 m s^{-1}'];
text(-0.5,-120,t,'fontweight','bold','fontsize',8)

print -depsc ../figuras/ex_metdinref.eps
!epstopdf ../figuras/ex_metdinref.eps









