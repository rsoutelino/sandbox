clear all;
clear all;close all;

N=2^10;
rt=2.312;

t=0:1:N-1;
t=t*rt;
w=2*pi/51/rt;
h=.25*pi;
A=1;
ys=A*sin(w*t+h);
yn=2*randn(1,N);
y=ys;%+yn;  

fy=fft(y); 

yr=zeros(1,N);

An = real(fy);
Bn = imag(fy);

Amp = sqrt(An.^2+Bn.^2);
Fase = atan2(Bn,An);

k=0:N-1;n=0:N-1;

stop

yr = (1/N)*sum(repmat(Amp,N,1)'.*cos(repmat(Fase,N,1)'+2*pi*(k'*n)/N));  

% a inversa (yr) é igual aos dados de entrada!
title('Aqui tem 2 curvas')
plot(t,y,t,yr,'y--')

% o vetor de amplitudes é 'espelhado' em relação ao zero
% o eixo das freqüências é construído levando-se isto em conta
fre=2*pi*(1:N/2)/t(N);
A=(2/N)*fliplr(Amp((N/2+1):N));
figure
subplot(2,1,1),plot(fre,A)
title('Amplitudes')
subplot(2,1,2),plot(fre,A);axis([0 .3 0 1]);
title('Zoom no Pico ')

pasbax1=double(fre<0.2);
filtro=[pasbax1 fliplr(pasbax1)];
Ampf=Amp.*real(filtro);
yrf1 = (1/N)*sum(repmat(Ampf,N,1)'.*cos(repmat(Fase,N,1)'+2*pi*(k'*n)/N));  

pasbax2=double(fre<0.1);
filtro=[pasbax2 fliplr(pasbax2)];
Ampf=Amp.*real(filtro);
yrf2 = (1/N)*sum(repmat(Ampf,N,1)'.*cos(repmat(Fase,N,1)'+2*pi*(k'*n)/N));  

pasban=double(fre>0.04 & fre<0.06);
filtro=[pasban fliplr(pasban)];
Ampf=Amp.*real(filtro);
yrf3 = (1/N)*sum(repmat(Ampf,N,1)'.*cos(repmat(Fase,N,1)'+2*pi*(k'*n)/N));  

figure
subplot(4,2,1),plot(t,y,'--b',t,ys,'r');axis([0 1000 -5 5])
title('Dados Originais e Sinal Senoidal')
subplot(4,2,4),plot(t,ys,'--b',t,yrf1,'r');axis([0 1000 -3 3])
title('Filtrado com Passa Baixa')
subplot(4,2,6),plot(t,ys,'--b',t,yrf2,'r');axis([0 1000 -2 2])
title('Filtrado com Passa Baixa')
subplot(4,2,8),plot(t,ys,'--b',t,yrf3,'r');axis([0 1000 -2 2])
title('Filtrado com Passa Banda')

subplot(4,2,3),plot(fre,A,fre,pasbax1.*A,fre,pasbax1);axis('tight')
title('Filtro Passa Baixa')
subplot(4,2,5),plot(fre,A,fre,pasbax2.*A,fre,pasbax2);axis('tight')
title('Filtro Passa Baixa')
subplot(4,2,7),plot(fre,A,fre,pasban.*A,fre,pasban);axis('tight')
title('Filtro Passa Banda')

