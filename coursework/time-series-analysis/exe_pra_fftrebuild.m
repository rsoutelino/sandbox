clear all;
clear all;close all;

N=2^10;

rt=2;
t=0:1:N-1;
t=t*rt;
w=2*pi/101/rt;
h=.25*pi;
A=1;
ys=A*sin(w*t+h);
yn=.3*randn(1,N);
y=yn+ys;  

fy=fft(y); 

% confira: a fft de y é complexa
whos fy

yc=zeros(1,N);
yr=zeros(1,N);
yr1=yr;yc2=yr;yc3=yr;

%Remova os comentários e veja como loop é ineficiente
%tic
%for n=1:N
%  for k=1:N
%    yc(k) = fy(k)*exp((2*pi*i)*(k-1)*(n-1)/N);
%  end
%  yr(n)=(1/N)*sum(yc);
%end
%toc

%tic
k=0:N-1;
for n=1:N

  % expressa como exponencial
  yc = fy.*exp((2*pi*i)*k*(n-1)/N);
  yr(n) = (1/N)*sum(yc);

  % expressa como produto de complexos
  yc1 = fy.*(cos(2*pi*k*(n-1)/N)+i*sin(2*pi*k*(n-1)/N));
  yr1(n) = (1/N)*sum(yc);
  
  % expressa como produto de reais
  An = real(fy);
  Bn = imag(fy);
  yc1r = An.*cos(2*pi*k*(n-1)/N);
  yc1i = Bn.*sin(2*pi*k*(n-1)/N);
  yc2(n) = (1/N)*sum(yc1r-yc1i);

  % note que o par An Bn contém informação sobre amplitude E fase
  % pois a fase é atan2(Bn,An) e a amplitude é sqrt(An^2+Bn^2)
  Amp = sqrt(An.^2+Bn.^2);
  Fase = atan2(Bn,An);
  yc3(n) = (1/N)*sum(Amp.*cos(Fase+2*pi*k*(n-1)/N));
  
  
  
end
%toc

subplot(2,1,1),plot(t,y,t,real(yr),t,real(yr1),t,real(yc2),t,real(yc3));
axis('tight')
subplot(2,1,2),plot(t,y,t,real(yr)+1,'.',t,real(yr1)+2,'-.',...
		    t,real(yc2)+3,':',t,real(yc3)+4,'--')
axis([888 1000 -2 5])



