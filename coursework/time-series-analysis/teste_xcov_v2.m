clear all;close all

load exe_pra_2.mat

t=t*rt;

y1=yns1;
y2=yns2;


% FFT EM CADA UMA DELAS
fy1=fft(y1);
fy2=fft(y2);

fy=fy2.*conj(fy1); 

fy=fliplr(fy((length(t)/2)+1:length(t)));

freq=2*pi*[1:length(t)/2]/(length(t)*rt);

xamp=sqrt(fy.*conj(fy));

xphase=atan2(imag(fy),real(fy));

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


subplot(4,1,1),loglog(freq,xamp);
hold on; 
% plot([w1 w1],[min(xamp) max(xamp)],'r') ; hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')

subplot(4,1,2),semilogx(freq,(180/pi)*xphase);
hold on;  
% plot([w1 w1],[-180 180],'r',[w1 w1 ],[-90 -90],'or') ;
hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Angulo')

subplot(4,1,3),loglog(freq,xamp);
hold on;  
%plot([w1 w1],[min(xamp) max(xamp)],'r') ; hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Potencia Espectral (amplitude)')

subplot(4,1,4),semilogx(freq,(180/pi)*xphase);
hold on; 
% plot([w1 w1],[-180 180],'r',[w1 w1 ],[-90 -90],'or') ; hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Angulo (fase)')




