clear all
close all

% O de sempre...
T=55; % Periodo da onda (ja sabemos de antemao)
t=0:1:1000; % tempo
N=length(t);

A0=1;
A1=3.75e-3;
A2=1.4;
A3=pi/4;
EA=0.75;

B0=5;
B1=-3e-3;
B2=1.5;
B3=-pi/4;
EB=0.65;

aa=A0+A1*t;
sa=A2*sin((2*pi/T)*t+A3);
na=EA*randn(1,N);
Ya=aa+sa+na;

bb=B0+B1*t;
sb=B2*sin((2*pi/T)*t+B3);
nb=EB*randn(1,N);
Yb=bb+sb+nb;
%Yb=aa+sa+nb;

subplot(2,2,1),plot(t,Ya,'b',t,Yb,'r')
title('The series')

Yam=mean(Ya);
Ybm=mean(Yb);
Yas=std(Ya);
Ybs=std(Yb);

% Estude E&T formula 5.3.3 and beyond

% C ---> Covariacia              R ---> coRrela��o 
%stop
for k=0:N-1

  Caa(k+1)=(1/(N-k))*sum( (Ya(1:N-k)-Yam).*(Ya(k+1:N)-Yam) );
  Cbb(k+1)=(1/(N-k))*sum( (Yb(1:N-k)-Ybm).*(Yb(k+1:N)-Ybm) );

  Raa(k+1)=(1/(N-k))*sum( Ya(1:N-k).*Ya(k+1:N) );
  Rbb(k+1)=(1/(N-k))*sum( Yb(1:N-k).*Yb(k+1:N) );

  Cab(k+1)=(1/(N-k))*sum( (Ya(1:N-k)-Yam).*(Yb(k+1:N)-Ybm) );
  Rab(k+1)=(1/(N-k))*sum( Ya(1:N-k).*Yb(k+1:N) );
    
end

tv=-fliplr(t);
tau=[tv(1:length(tv)-1) t];

aaC=fliplr(Caa);
aaR=fliplr(Raa);
Caa=[aaC(1:length(aaC)-1) Caa];
Raa=[aaR(1:length(aaR)-1) Raa];

bbC=fliplr(Cbb);
bbR=fliplr(Rbb);
Cbb=[bbC(1:length(bbC)-1) Cbb];
Rbb=[bbR(1:length(bbR)-1) Rbb];

baC=fliplr(Cab);
baR=fliplr(Rab);
Cab=[baC(1:length(baC)-1) Cab];
Rab=[baR(1:length(baR)-1) Rab];


subplot(2,2,2),plot(tau,Caa,'b',tau,Cbb,'r')
hold on
plot(0,Yas^2,'ob',0,Ybs^2,'or')
hold off
title('As Covariancias')

subplot(2,2,3),plot(tau,Raa,'b',tau,Rbb,'r')
hold on
plot(0,Yas^2+Yam^2,'+b',0,Ybs^2+Ybm^2,'+r')
hold off
title('As Autocorrelacoes')

subplot(2,2,4),plot(tau,Cab,'b',tau,Rab,'r')
title('As Covariancias e Correlacoes Cruzadas')

h =findobj('color',[1 0 0])
set(h,'color',[0.7 0.7 0.7])
h =findobj('color',[0 0 1])
set(h,'color',[.15 .15 .15])

% Exerc�cio: plote as fun��es NORMALIZADAS cf. pg. 374-376h =findobj('color',[0 0 1])
