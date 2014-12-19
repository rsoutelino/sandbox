%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Rafael Guarino Soutelino %%%%%%%%%%%%%%%%%
%%%% Lista 1, problema 4 - Oceanografia Observacional %%%
%%%%%  		     Maio / 2006		     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;close all;

load exe_pra_2.mat


t=t*rt;

npf=31; % tamanho do filtro

inxs=(npf-1)/2;

% FILTRO
myfilter=blackman(npf);

myfilter=myfilter/sum(myfilter);

nrep=20000; %pedacos da serie temporal

%%% DIVIDE A SERIE EM PEDACOS E TIRA A MEDIA

k=0;
for j=nrep:nrep:length(t)
    i=j-nrep+1;
    k=k+1;
    
  y1=conv(yns1(i:j),myfilter);
  y2=conv(yns2(i:j),myfilter);
  
  y1=y1(inxs+1:length(y1)-inxs);
  y2=y2(inxs+1:length(y2)-inxs);
  
  fy1=fft(y1); 
  fy2=fft(y2); 

  fy=fy2.*conj(fy1);
  fy=fliplr(fy((length(t(i:j))/2)+1:length(t(i:j))));
  
  freq=2*pi*[1:length(t(i:j))/2]/(length(t(i:j))*rt);

  xamp(k,:)=sqrt(fy.*conj(fy));

  xphase(k,:)=atan2(imag(fy),real(fy));
  
  if j==nrep,
    sfy=fy;
  else
    sfy=sfy+fy;
  end
end

mamp=mean(xamp);
mphase=mean(xphase);

dpamp=std(xamp)./(sqrt(k));
dpphase=std(xphase)./(sqrt(k));

dpampup=mamp+dpamp;
dpampdown=mamp-dpamp;

% dp1=mfy+dp;
% dp2=mfy-dp;

Pico=near(mamp,max(mamp),2);

Freq=freq(Pico);

% PLOTAGEM

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','A4',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');


subplot(3,1,1),

plot(t(1:500),yns1(1:500),'r');
hold on;
plot(t(1:500),yns2(1:500),'k');
axis('tight');
title('Series brutas')
xlabel('tempo')
legend('Serie 1','Serie 2')

subplot(3,1,2), 

loglog(freq,dpampup,'c');grid on; hold on
loglog(freq,dpampdown,'c');
loglog(freq,mamp,'k');
axis('tight');% axis([2e-3 2 10e0 10e7]);
xlabel('Frequencia')
ylabel('Potencia Espectral')
legend('Espectro cruzado','Desvio Padrao')

subplot(3,1,3),

semilogx(freq,(180/pi)*mphase,'k');
axis('tight'); 
grid on;
xlabel('Frequencia')
ylabel('Angulo')
legend('Fase')

print -depsc prob_4_1.eps
!epstopdf prob_4_1.eps

figure(2) 
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','A4',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');


subplot(4,1,1),

loglog(freq,mamp,'k');
axis('tight'); % axis([0.10 0.35 25e5 .5e8]);
grid on;
hold on; 
plot([Freq(1) Freq(1)],[min(mamp) max(mamp)],'k');
title('Zoom pico 1')
ylabel('Amplitude')

subplot(4,1,2),

semilogx(freq,(180/pi)*mphase);
axis('tight'); %axis([0.10 0.35 -180 180]);
grid on; hold on;
plot([Freq(1) Freq(1)],[-180 180],'k',[Freq(1) Freq(1)],[-90 -90],'ok') ;
ylabel('Fase')

subplot(4,1,3),

loglog(freq,mamp,'k');
axis('tight'); %axis([0.04 0.15 30e5 .5e8]);
grid on; hold on;
plot([Freq(2) Freq(2)],[min(mamp) max(mamp)],'k');
title('Zoom pico 2')
ylabel('Amplitude')

subplot(4,1,4),
semilogx(freq,(180/pi)*mphase);
axis('tight'); %axis([0.04 0.15 -180 180]);
grid on; hold on;
plot([Freq(2) Freq(2)],[-180 180],'k',[Freq(2) Freq(2)],[90 90],'ok') ;
xlabel('Frequencia')
ylabel('Fase')

print -depsc prob_4_2.eps
!epstopdf prob_4_2.eps