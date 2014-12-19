clear all;close all;clc

g = 10;
H = 100;
lambda = 100;
k = 2*pi/lambda;
w = sqrt(g*k*tanh(k*H));

x = 1:2:200;

%  eta1 = cos(k*x - w*t);
%  eta2 = cos(-k*x - w*t);

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 9],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 9],...
        'ShareColors','off',...
        'Clipping','on');

subplot(611)
title('Onda Estacionaria','fontweight','bold')
hold on

for t = 0:5
   if t==5,break,end
   eta = 2*cos(k*x)*cos(w*t);
   eval(['subplot(6,1,',num2str(t+1),')'])
   plot(x,eta,'linewidth',2); hold on
   plot(x,zeros(size(x)),'k--')
   axis([x(1) x(end) -2.5 2.5])
   set(gca,'yticklabel',[])
   set(gca,'xticklabel',[])
   tit = ['t = ',num2str(t)];
   ylabel(tit,'fontweight','bold')
end

% representando vetores de velocidade em t=0

t=1; 
eta = 2*cos(k*x)*cos(w*t);
util = ((k*g)/(w*cosh(k*H)))*cosh(k*H)*(cos(k*x - w*t)-cos(-k*x - w*t));
wtil = ((2*k*g)/(w*cosh(k*H)))*sinh(k*H)*(sin(k*x - w*t)-sin(-k*x - w*t));

subplot(616)
plot(x,zeros(size(x)),'k--'); hold on
quiver(x,eta,util,wtil,0.05,'k');hold on
plot(x,eta,'r','linewidth',2)
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])
ylabel('t=0','fontweight','bold')


print -dpng exe5.png
%  !epstopdf exe5.eps
%  !rm -rf exe5.eps







