%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%    Lista Pratica # 7 - Problema  1                                       %
%    Prof. Ilson Carlos Almeida da Silveira                                %
%                                                                          %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all ; clear all ; clc

%load saidas do programa da lista 5
load F.mat;
load medios.mat

Rd_0=15583.3942*1e3;
Rd_1=196.132203*1e3;
Rd_2=124.124119*1e3;
Rd_3=76.5253272*1e3;
Rd_4=56.8374664*1e3;
Rd_5=44.8154785*1e3;

g=9.8;
H=4000;
f0=2*7.3e-5*sin(5*pi/180);

%figure
%plot(dens(1:4000), -pm(1:4000), 'linewidth', 2);
%hold on
%plot(dens_est, -pm(1:4000), '--r', 'linewidth', 2)
%set(gca,'fontsize',16)
%title('Pefil de Densidade','fontsize',20)
%ylabel('Profundidade (m)','fontsize',18)
%xlabel('Dendidade (kg/m3)', 'fontsize',18)
%axis('square')

%---------------------------------------
%Calibracao(1): Bombeamento de Ekman
%---------------------------------------

hb1=H/(1+F1(1)^2);
epsib=(Rd_1^2*f0^2*H)/(g*hb1*(H-hb1));

%---------------------------------------
%Calibracao (2):Interacoes nao-lineares
%---------------------------------------

qsi=(1/H)*10*sum(F1.^3);
gama=(1/4)*((sqrt(qsi^2+4)-qsi))^2;
hn1=(gama*H)/(1+gama);
epsin=(Rd_1^2*f0^2*H)/(g*hn1*(H-hn1));

hb=[hb1 H-hb1]; HB=round(hb1/10);
hn=[hn1 H-hn1]; HN=round(hn1/10);
fb=[sqrt((H-hb1)/hb1); -sqrt(hb1/(H-hb1))];
fn=[sqrt((H-hn1)/hn1); -sqrt(hn1/(H-hn1))];

% calcula densidade em cada camada

D0=mean(Dm);
Dm=Dm-1000;

rho2b=mean(Dm(4000-HB:4000));
rho1b=rho2b-D0*epsib;

rho2n=mean(Dm(4000-HN:4000));
rho1n=rho2n-D0*epsin;

% constroi perfis dos modos discretos e densidade

Fn=zeros(400,1);
Fb=Fn;Rhob=Fn;Rhon=Fn;

Fb(1:HB)=ones(HB,1)*fb(1);
Fb(HB+1:400)=ones(400-HB,1)*fb(2);

Rhob(1:HB)=ones(HB,1)*rho1b;
Rhob(HB+1:400)=ones(400-HB,1)*rho2b;

Fn(1:HN)=ones(HN,1)*fn(1);
Fn(HN+1:400)=ones(400-HN,1)*fn(2);

Rhon(1:HN)=ones(HN,1)*rho1n;
Rhon(HN+1:400)=ones(400-HN,1)*rho2n;


figure(1)
     subplot(1,2,1)
plot(Rhob,-(0:10:3990),'--c',Dm,-P,'m','linewidth',3)
     legend('2-camadas','WESTRAX 4')
title('Densidade Modelada')
xlabel('Densidade [kg m^3]')
ylabel('Profundidade (m)')


     subplot(1,2,2)
plot([0 0],[0 -3990],'k')
hold on
plot(Fb,-(0:10:3990),'--b',[1 1],[0 -3990],'--g','linewidth',3)
     axis([-1 4 -3990 0])
title('Modos Normais - Bombeamento de Ekman')
xlabel('Amplitude')
ylabel('Profundidade (m)')

% print -depsc2 modos_b.eps
 
 figure(2)
     subplot(1,2,1)
plot(Rhon,-(0:10:3990),'--c',Dm,-P,'m','linewidth',3)
     legend('2-camadas','WESTRAX 4')
title('Densidade Modelada')
xlabel('Densidade [kg m^3]')
ylabel('Profundidade (m)')

     subplot(1,2,2)
 plot(Fn,-(0:10:3990),'--r',[1 1],[0 -3990],'--g','linewidth',3)
 hold on
 plot([0 0],[0 -3990],'k')
     axis([-1 4 -3990 0])
 title('Modos Normais - Intera�es N� Lineares')
 xlabel('Amplitude')
 ylabel('Profundidade (m)')
 set(gcf,'PaperPositionMode','auto')

% print -depsc2 modos_i.eps

save param_naolinear epsin hn
 
 


