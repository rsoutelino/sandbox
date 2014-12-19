clear all;close all;

rt=1;

t=0:1:2^(16)-1;
t=t*rt;


w=2*pi/628;
h=.25*pi;
A=0.5;
ys=A*sin(w*t+h);

freq=2*pi*[1:length(t)/2]/max(t);

npf=331;

inxs=(npf-1)/2;

myfilter=ones(1,npf);

%myfilter=triang(npf);

%myfilter=hamming(npf);

%myfilter=blackman(npf);


myfilter=myfilter/sum(myfilter);

nrep=30;

for j=1:nrep,
  yn=10*A*randn(size(t));
  
  y=conv(yn+ys,myfilter);
  
  y=y(inxs+1:length(y)-inxs);
  
  fy=fft(y); fy1=fy;
  fy=fy.*conj(fy); fy2=fy;
  fy=fliplr(fy((length(t)/2)+1:length(t)));
  if j==1,
    sfy=fy;
  else
    sfy=sfy+fy;
  end
end

mfy=sfy/nrep;



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

subplot(3,1,1),plot(t(1:2000),ys(1:2000));axis('tight')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais')

subplot(3,1,2),plot(t(1:2000),y(1:2000),'k',t(1:2000),ys(1:2000),'c');axis('tight')
legend('y','ys')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais + Ruido')

subplot(3,1,3),loglog(freq,mfy);axis([6e-4 .1 1 1e8]);grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')



