clear all;close all;clc

% criando vetor de tempo
t = 0:0.1:50;
w = 0.44;
H = 100;

% plotando eta
eta = 0.5*cos(w*t);

figure(1)
set(gcf,'color','w')
plot(t,zeros(size(t)),'--k'); hold on
plot(t,eta,'linewidth',2)
axis([0 t(end) -4 1])
xlabel('t [s]','fontweight','bold')
title('Elevacao da Superficie Livre do Mar','fontweight','bold')
ylabel('\eta [m]','fontweight','bold')
pbaspect([1 0.4 1])
print -depsc eta.eps
!epstopdf eta.eps
!rm -rf eta.eps

% plotando ptil nos niveis z
p0 = 5127*cos(w*t);
p10 = 4225*cos(w*t);
p50 = 2099*cos(w*t);
p100 = 1363*cos(w*t);

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
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(4,1,1)
plot(t,p0,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -10000 10000])
title('Perturbacao de pressao','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = 0 m')

subplot(4,1,2)
plot(t,p10,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -10000 10000])
set(gca,'xticklabel',[])
legend('z = -10 m')

subplot(4,1,3)
plot(t,p50,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -10000 10000])
ylabel('perturbacao de pressao [N m^{-2}]','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = -50 m')

subplot(4,1,4)
plot(t,p100,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -10000 10000])
legend('z = -100 m')
xlabel('t [s]','fontweight','bold')

print -depsc p_til.eps
!epstopdf p_til.eps
!rm -rf p_til.eps

% plotando util nos niveis z

u0 = 0.22*cos(w*t);
u10 = 0.18*cos(w*t);
u50 = 0.09*cos(w*t);
u100 = 0.06*cos(w*t);

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
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(4,1,1)
plot(t,u0,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
title('Perturbacao de velocidade horizontal [u]','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = 0 m')

subplot(4,1,2)
plot(t,u10,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
set(gca,'xticklabel',[])
legend('z = -10 m')

subplot(4,1,3)
plot(t,u50,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
ylabel('perturbacao de velocidade [m s^{-1}]','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = -50 m')

subplot(4,1,4)
plot(t,u100,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
legend('z = -100 m')
xlabel('t [s]','fontweight','bold')

print -depsc u_til.eps
!epstopdf u_til.eps
!rm -rf u_til.eps


% plotando wtil nos niveis z

w0 = 0.22*sin(w*t);
w10 = 0.17*sin(w*t);
w50 = 0.07*sin(w*t);
w100 = 0*sin(w*t);


figure(4)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(4,1,1)
plot(t,w0,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
title('Perturbacao de velocidade vertical [w]','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = 0 m')

subplot(4,1,2)
plot(t,w10,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
set(gca,'xticklabel',[])
legend('z = -10 m')

subplot(4,1,3)
plot(t,w50,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
ylabel('perturbacao de velocidade [m s^{-1}]','fontweight','bold')
set(gca,'xticklabel',[])
legend('z = -50 m')

subplot(4,1,4)
plot(t,w100,'linewidth',2); hold on
plot(t,zeros(size(t)),'--k'); hold on
axis([0 t(end) -0.5 0.5])
legend('z = -100 m')
xlabel('t [s]','fontweight','bold')

print -depsc w_til.eps
!epstopdf w_til.eps
!rm -rf w_til.eps





















