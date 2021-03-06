clear all;close all;

rt=1; % resolucao temporal

t=0:rt:10000-1;

w1=2*pi/2344; % regula a frequencia
h1=.25*pi; % regula a fase
A1=0.8; % regula a amplitude
ys1=A1*sin(w1*t+h1); % montando o sinal

w2=2*pi/718;
h2=.75*pi;
A2=0.6;
ys2=A2*sin(w2*t+h2);

w3=2*pi/309;
h3=.25*pi;
A3=0.5;
ys3=A3*sin(w3*t+h3);

ys=ys1+ys2+ys3; % soma dos tres sinais

A=2*sqrt(A1^2+A2^2+A3^2); % criando uma amplitude para o sinal aleatorio

yn=A*randn(size(t)); % criando um sinal aleatorio com distribuicao normal

y=yn+ys; % somando o sinal aleatorio aos tres sinais previamente criados 

fy=fft(y); fy1=fy; % transformada de fourrier

fy=fy.*conj(fy); fy2=fy; % eliminando a parte imaginaria

fy=fliplr(fy((length(t)/2)+1:length(t))); % tirando a metade e invertendo

freq=2*pi*[1:length(t)/2]/max(t); % montando eixo de frequencias
f1=[-fliplr(freq) freq];

pow1=(std(y)*10000).^2;
pow2=2*sum(fy);
[pow1 pow2]

figure(1)
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

subplot(4,2,7),loglog(freq,fy);axis('tight');grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')

subplot(4,2,8),loglog(freq,fy);axis([6e-4 .05 1000 2e7]);grid on
hold on;plot([w1 w2 w3],[1e6 1e6 1e6],'or');
plot([w1 w2 w3],[1e6 1e6 1e6],'+k');hold off
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')

%%% APLICANDO MÉDIA MOVEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yf=weim(1001,'hann',y);


fyf=fft(yf); fy1=fyf;

fyf=fyf.*conj(fyf); fy2f=fyf;

fyf=fliplr(fyf((length(t)/2)+1:length(t)));

freq=2*pi*[1:length(t)/2]/max(t);
f1=[-fliplr(freq) freq];

pow1=(std(yf)*10000).^2;
pow2=2*sum(fyf);
[pow1 pow2]

figure(2)
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

subplot(4,1,1),plot(t,ys1,'-',t,ys2,'--',t,ys3,'-.')
legend('ys1','ys2','ys3')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais')

subplot(4,1,2),plot(t,yf,'k',t,ys,'c')
legend('y','ys')
xlabel('Tempo')
ylabel('Velocidade')
title('Sinais + Ruido')

subplot(4,2,5),plot(f1,real(fy1),'b',f1,imag(fy1),'r');grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('Partes real e imaginaria da fft de y')
legend('real','imag')

subplot(4,2,6),semilogy(f1,fy2f);grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('Modulo da fft de y')

subplot(4,2,7),loglog(freq,fyf);axis('tight');grid on
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')

subplot(4,2,8),loglog(freq,fyf);axis([6e-4 .05 1000 2e7]);grid on
hold on;plot([w1 w2 w3],[1e6 1e6 1e6],'or');
plot([w1 w2 w3],[1e6 1e6 1e6],'+k');hold off
xlabel('Frequencia')
ylabel('Potencia Espectral')
title('PSD')











