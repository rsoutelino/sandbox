clear all;close all;format long g

load exe_pra_2.mat;

t=t*rt;

%%% CONDUZINDO ANALISE ESPECTRAL DAS SERIES 
npf=61; % tamanho do filtro

inxs=(npf-1)/2;

% analise espectral das series completas brutas

fy1=fft(yns1); 
fy1=fy1.*conj(fy1); 
fy1=fliplr(fy1((length(t)/2)+1:length(t)));

fy2=fft(yns2); 
fy2=fy2.*conj(fy2); 
fy2=fliplr(fy2((length(t)/2)+1:length(t)));



%%% APLICACAO DO FILTRO %%%%%%%%%%%%%%%%%%

inxs=(npf-1)/2;

% myfilter=ones(1,npf);
% myfilter=triang(npf);
% myfilter=hamming(npf); % tipo de filtro 
myfilter=blackman(npf);

myfilter=myfilter/sum(myfilter);
y1=conv(yns1,myfilter);
y2=conv(yns2,myfilter);
y1=y1(inxs+1:length(y1)-inxs);
y2=y2(inxs+1:length(y2)-inxs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% analise espectral da serie completa filtrada
   
  ffy1=fft(y1); 
  ffy1=ffy1.*conj(ffy1); 
  ffy1=fliplr(ffy1((length(t)/2)+1:length(t)));

  ffy2=fft(y2); 
  ffy2=ffy2.*conj(ffy2); 
  ffy2=fliplr(ffy2((length(t)/2)+1:length(t)));

dp1=[];
dp2=[]; % para guardar o std em cada segmento do espectro

nrep=40000; % numero de amostras em cada segmento do espectro

freq1=2*pi*[1:length(t)/2]/max(t); % eixo frequencia -espectro completo
freq=2*pi*[1:length(t(1:nrep))/2]/max(t(1:nrep)); % espectro medio


% LOOP PARA CALCULAR O ESPECTRO EM CADA SEGMENTO E FAZER A MEDIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=0; % iniciando contador
for j=nrep:nrep:length(t)
  n=j-nrep+1;
  k=k+1;

  my1=conv(yns1(n:j),myfilter);
  my1=my1(inxs+1:length(my1)-inxs);
   
  my1=fft(my1);
  my1=my1.*conj(my1);
  my1=fliplr(my1((length(t(n:j))/2)+1:length(t(n:j))));

  my2=conv(yns2(n:j),myfilter);
  my2=my2(inxs+1:length(my2)-inxs);
   
  my2=fft(my2);
  my2=my2.*conj(my2);
  my2=fliplr(my2((length(t(n:j))/2)+1:length(t(n:j))));

  if j==nrep,
    sfy1=my1;
    sfy2=my2;
  else
    sfy1=sfy1+my1;
    sfy2=sfy2+my2;
  end
  dp1(k,:)=my1;
  dp2(k,:)=my2; % matriz qu % matriz que guarda todos os espectros para posterior calculo do desvio padrao
end

mfy1=sfy1/nrep;
mfy2=sfy2/nrep;

dpd1=std(dp1); % dpd=dpd./(sqrt(length(t)/nrep));
dpd2=std(dp2);

%%% PLOTAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','landscape',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 15 11],...
        'ShareColors','off',...
        'Clipping','on');


% serie bruta + filtrada
subplot(4,2,1),plot(t(1:nrep),yns1(1:nrep),'r',t(1:nrep),y1(1:nrep),'k')
axis('tight')
legend('Serie 1 Bruta','Serie Filtrada')
title('Serie 1','fontsize',12,'fontweight','bold')

subplot(4,2,2),plot(t(1:nrep),yns2(1:nrep),'r',t(1:nrep),y2(1:nrep),'k')
axis('tight')
legend('Serie 2 Bruta','Serie Filtrada')
title('Serie 2','fontsize',12,'fontweight','bold')

%------------------------------------------------------------------------

% espectro bruto
subplot(4,2,3),loglog(freq1,fy1,'r');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Bruto - Serie 1')

subplot(4,2,4),loglog(freq1,fy2,'r');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Bruto - Serie 2')

%------------------------------------------------------------------------

% espectro filtrado
subplot(4,2,5),loglog(freq1,ffy1,'k');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Filtrado - Serie 1')

subplot(4,2,6),loglog(freq1,ffy2,'k');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Filtrado - Serie 2')

%------------------------------------------------------------------------

% espectro medio filtrado
subplot(4,2,7),loglog(freq,mfy1);axis([1e-4 1e-1 1e-2 1e6]);grid on
hold on
%  loglog(freq,mfy1+dpd1,'c')
%  loglog(freq,mfy1-dpd1,'k')
%  axis('tight')
xlabel('Frequencia')
ylabel('Potencia Espectral')
title(['Espectro Medio da serie 1 com ',num2str(floor(length(t)/nrep)),' segmentos'])

subplot(4,2,8),loglog(freq,mfy2);axis([1e-4 1e-1 1e-2 1e6]);grid on
hold on
%  loglog(freq,mfy2+dpd2,'c')
%  loglog(freq,mfy2-dpd2,'k')
%  axis('tight')
xlabel('Frequencia')
ylabel('Potencia Espectral')
title(['Espectro Medio da serie 2 com ',num2str(floor(length(t)/nrep)),' segmentos'])


%------------------------------------------------------------------------
%  
%  figure(2)
%  
%  plot(t,yns1,'r',t,yns2,'k');axis([1 3000 -400 400])

% SALVANDO AS FIGURAS
print -depsc prob_4.eps
!epstopdf prob_4.eps

%%% CHECANDO COEF DE CORRELACAO E ERRO MEDIO QUADRATICO

[r,p,rlo,rup]=corrcoef(yns1, yns2);

corr=r(1,2)
prob=p(1,2)

rms= sqrt(sum( (yns1-yns2).^2)./sum(yns1.^2))

