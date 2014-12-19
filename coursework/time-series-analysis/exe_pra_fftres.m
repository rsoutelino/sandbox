
clear all;close all;
m=round(logspace(2,5,10));


for n=1:length(m)
N=m(n)+mod(m(n),2);
  %N=1e5;
rt=3.6;
A=[1 2 3 1];
T=[1900 1400 1100]; w=2*pi./T;
h=pi./[2 3 4];

t=0:1:N-1;
t=t*rt;

ys1=A(1)*sin(w(1)*t+h(1));
ys2=A(2)*sin(w(2)*t+h(2));
ys3=A(3)*sin(w(3)*t+h(3));
yn=A(4)*randn(1,N);

ys=ys1+ys2+ys3;
y=yn+ys;

fy=fft(y);

Am=fy.*conj(fy);
Am=2*fliplr(Am((N/2+1):N))/(N^2);
freq=2*pi*[1:N/2]/t(N);
f1=[-fliplr(freq) freq];

loglog(freq,Am,'-b',freq,Am,'pr',w,[std(ys1)^2 std(ys2)^2 std(ys3)^2],'ok');
grid on
set(gca,'xlim', [min(w)/3 max(w)*3],'ylim',[1e-5 10],'xtick',...
	logspace(-5,1))
xlabel('Frequencia')
ylabel('Potencia Espectral')
title(['PSD baseada em ' num2str(N) ' pontos'])
drawnow;pause
end
