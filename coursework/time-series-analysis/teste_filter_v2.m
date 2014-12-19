%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  Rafael Guarino Soutelino %%%%%%%%%%%%%%%%%
%%%% Lista 1, problema 3 - Oceanografia Observacional %%%
%%%%%  		     Maio / 2006		     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

load exe_pra_1.mat;

%  load yns.dat % serie sintetica, comentar para nao usar

t=t*rt;

npf=61; % TAMANHO DO FILTRO

inxs=(npf-1)/2;

% analise espectral da serie completa bruta

fy0=fft(yns); 
fy0=fy0.*conj(fy0); 
fy0=fliplr(fy0((length(t)/2)+1:length(t)));

%%% FILTROS %%%%%%%%%%%%%%%%%%

%  myfilter=ones(1,npf);
%  myfilter=triang(npf);
%  myfilter=hamming(npf);
myfilter=blackman(npf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myfilter=myfilter/sum(myfilter);

%%% analise espectral da serie completa filtrada

y1=conv(yns,myfilter);

  y1=y1(inxs+1:length(y1)-inxs);
   
  fy_1=fft(y1); fy_2=fy_1;
  fy_1=fy_1.*conj(fy_1); fy_3=fy_1;
  fy_1=fliplr(fy_1((length(t)/2)+1:length(t)));


dpd=[]; % para guardar o std em cada segmento do espectro
nrep=20000; % numero de amostras em cada segmento do espectro

freq1=2*pi*[1:length(t)/2]/max(t); % eixo frequencia -espectro completo
freq=2*pi*[1:length(t(1:nrep))/2]/max(t(1:nrep)); % espectro medio


% LOOP PARA CALCULAR O ESPECTRO EM CADA SEGMENTO E FAZER A MEDIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=0; % iniciando contador
for j=nrep:nrep:length(t)
  n=j-nrep+1;
  k=k+1;

  y=conv(yns(n:j),myfilter);

  y=y(inxs+1:length(y)-inxs);
   
  fy=fft(y); fy1=fy;
  fy=fy.*conj(fy); fy2=fy;
  fy=fliplr(fy((length(t(n:j))/2)+1:length(t(n:j))));

  dp(k,:)=fy; % matriz que guarda todos os espectros para posterior calculo do espectro medio e seu desvio padrao em cada frequencia
end

mfy=mean(dp);

dpd=std(dp);
err=mfy+dpd;

%%% PLOTAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','A4',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');


% serie bruta + filtrada
subplot(4,1,1),plot(t(1:nrep),yns(1:nrep),'r',t(1:nrep),y1(1:nrep),'k')
axis('tight')
legend('Serie Bruta','Serie Filtrada')
title('Analise Espectral','fontsize',12,'fontweight','bold')

% espectro bruto
subplot(4,1,2),loglog(freq1,fy0,'r');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Bruto')

% espectro filtrado
subplot(4,1,3),loglog(freq1,fy_1,'b');axis([1e-4 1e-1 1e-2 1e12]);grid on
ylabel('Potencia Espectral')
title('Espectro Total Filtrado')

% espectro medio filtrado + desvio padrao
subplot(4,1,4),loglog(freq,err,'c');axis([1e-4 1e-1 1 1e10]);grid on
hold on
loglog(freq,mfy,'k')
legend('Espectro + dp','Espectro')
xlabel('Frequencia')
ylabel('Potencia Espectral')
title(['Espectro Medio com ',num2str(floor(length(t)/nrep)),' segmentos'])


% SALVANDO AS FIGURAS
print -depsc prob_3.eps
!epstopdf prob_3.eps











