clear all;close all;

load exe_pra_2.mat


t=t*rt;

npf=61; % tamanho do filtro

inxs=(npf-1)/2;

% FILTRO
myfilter=hamming(npf);

myfilter=myfilter/sum(myfilter);

nrep=20000; %pedacos da serie temporal

% DIVIDE A SERIE EM PEDACOS E TIRA A MEDIA

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
subplot(3,1,1),plot(t(1:500),yns1(1:500),'r');hold on;plot(t(1:500),yns2(1:500),'k');axis('tight');
title('Series brutas')
xlabel('tempo')
legend('serie 1','serie 2')

subplot(3,1,2),loglog(freq,mamp,'k');axis('tight');axis([2e-3 2 10e0 10e7]); grid on;
hold on;
loglog(freq,dpampup,'r');
loglog(freq,dpampdown,'r');
xlabel('Frequencia')
ylabel('Potencia Espectral')
legend('espectro cruzado','desvio padrao')

subplot(3,1,3),semilogx(freq,(180/pi)*mphase,'r');axis('tight'); grid on;
xlabel('Frequencia')
ylabel('Angulo')
legend('fase')

figure(2)
subplot(4,1,1),loglog(freq,mamp,'r');axis('tight');axis([0.10 0.35 25e5 .5e8]);grid on;
hold on; plot([Freq(1) Freq(1)],[min(mamp) max(mamp)],'k');
title('Zoom pico 1')
ylabel('Potencia Espectral')

subplot(4,1,2),semilogx(freq,(180/pi)*mphase);axis('tight');axis([0.10 0.35 -180 180]);
grid on; hold on;plot([Freq(1) Freq(1)],[-180 180],'k',[Freq(1) Freq(1)],[-90 -90],'ok') ;
ylabel('Angulo')

subplot(4,1,3),loglog(freq,mamp,'r');axis('tight');axis([0.04 0.15 30e5 .5e8]);grid on;
hold on;plot([Freq(2) Freq(2)],[min(mamp) max(mamp)],'k');
title('Zoom pico 2')
ylabel('Potencia Espectral')

subplot(4,1,4),semilogx(freq,(180/pi)*mphase);axis('tight');axis([0.04 0.15 -180 180]);
grid on; hold on;plot([Freq(2) Freq(2)],[-180 180],'k',[Freq(2) Freq(2)],[90 90],'ok') ;
xlabel('Frequencia')
ylabel('Angulo')
