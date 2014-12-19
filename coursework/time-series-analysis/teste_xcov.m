clear all;clf

rt=.32;

t=0:1:2^14-1;
t=t*rt;


w1=2*pi/75;
h1=.25*pi;
A1=0.8;
An1=2.2;
y1=A1*sin(w1*t+h1)+An1*randn(size(t));

w2=w1;
h2=.75*pi;
A2=0.8;
An2=3.2;
y2=A2*sin(w2*t+h2)+An2*randn(size(t));

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
hold on;  plot([w1 w1],[min(xamp) max(xamp)],'r') ; hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')

subplot(4,1,2),semilogx(freq,(180/pi)*xphase);
hold on;  plot([w1 w1],[-180 180],'r',[w1 w1 ],[-90 -90],'or') ; hold off
axis('tight')
grid on
xlabel('Frequencia')
ylabel('Angulo')

subplot(4,1,3),loglog(freq,xamp);
hold on;  plot([w1 w1],[min(xamp) max(xamp)],'r') ; hold off
axis([5e-2 .2 1e3 3e7])
grid on
xlabel('Frequencia')
ylabel('Potencia Espectral (amplitude)')

subplot(4,1,4),semilogx(freq,(180/pi)*xphase);
hold on;  plot([w1 w1],[-180 180],'r',[w1 w1 ],[-90 -90],'or') ; hold off
axis([5e-2 .2 -180 180])
grid on
xlabel('Frequencia')
ylabel('Angulo (fase)')




