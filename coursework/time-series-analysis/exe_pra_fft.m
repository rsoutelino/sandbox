clear all;close all;

rt=1;

t=0:1:10000-1;
t=t*rt;


w1=2*pi/2321;
h1=.25*pi;
A1=0.8;
ys1=A1*sin(w1*t+h1);

w2=2*pi/718;
h2=.75*pi;
A2=0.6;
ys2=A2*sin(w2*t+h2);

w3=2*pi/309;
h3=.25*pi;
A3=0.5;
ys3=A3*sin(w3*t+h3);

ys=ys1+ys2+ys3;

A=2*sqrt(A1^2+A2^2+A3^2);

yn=A*randn(size(t));

y=yn+ys;

fy=fft(y); fy1=fy;

fy=fy.*conj(fy); fy2=fy;

fy=fliplr(fy((length(t)/2)+1:length(t)));

freq=2*pi*[1:length(t)/2]/max(t);
f1=[-fliplr(freq) freq];

pow1=(std(y)*10000).^2;
pow2=2*sum(fy);
[pow1 pow2]


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
        'Clipping','on');

subplot(4,1,1),plot(t,ys1,'-',t,ys2,'--',t,ys3,'-.')
legend('ys1','ys2','ys3')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais')

subplot(4,1,2),plot(t,y,'k',t,ys,'c')
legend('y','ys')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais + Ruido')

subplot(4,2,5),plot(f1,real(fy1),'b',f1,imag(fy1),'r');grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('Partes real e imaginaria da fft de y')
legend('real','imag')

subplot(4,2,6),semilogy(f1,fy2);grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('Modulo da fft de y')

subplot(4,2,7),loglog(freq,2*fy/(length(t)^2),...
		      [w1 w2 w3],[std(ys1)^2 std(ys2)^2 std(ys3)^2],'o');
axis('tight');grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')

subplot(4,2,8),loglog(freq,2*fy/(length(t)^2),...
		      [w1 w2 w3],[std(ys1)^2 std(ys2)^2 std(ys3)^2],'o');
axis([6e-4 .05 5e-5 .5]);grid on
set(gca,'ytick',[1e-3 1e-2 1e-1 1],'xtick',[1e-3 1e-2]) 
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')

ha=findobj(gcf,'type','axes');
set(ha,'fontsize',14)
ht=cell2mat(get(ha,'title'));
set(ht,'fontsize',16)
ht=cell2mat(get(ha,'xlabel'));
set(ht,'fontsize',14)
ht=cell2mat(get(ha,'ylabel'));
set(ht,'fontsize',14)
