%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Analise das series temporais - Sao sebastiao %%%%%%%%%%%%
%%%%%%%%%              Junho - 2006                 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% espectros de vento, corrente e mare
% correlacao entre as variaveis
% espectros cruzados


clear all;close all;clc
format short g

%  carrega as matrizes de cada serie, provenientes dos programas de leitura
load corr.mat
load mare.mat
load estmet.mat

%  Corrente Residual
CRu=corr_u-pout_corr_u';
CRv=corr_v-pout_corr_v';

%  Elevacao Residual
ElR=mare_amp-pout_mare';

%  Vetor de DATA Completo
DATAC=datenum('17-Feb-2006 19:00'):datenum('00:30:00'):datenum('23-Mar-2006 10:30:00');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          PLOTAGEM INICIAL              %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; set(gcf,'color','w')
subplot(511)
plot(DATAC,corr_u)
title('Corrente Zonal [m.s^{-1}]')
datetick('x',19)

subplot(512)
plot(DATAC,corr_v)
title('Corrente Meridional [m.s^{-1}]')
datetick('x',19)

subplot(513)
plot(DATAC,estmet_u)
title('Vento Zonal [m.s^{-1}]')
datetick('x',19)

subplot(514)
plot(DATAC,estmet_v)
title('Vento Meridional [m.s^{-1}]')
datetick('x',19)

subplot(515)
plot(DATAC,mare_amp)
title('Amplitude de Mare [m]')
datetick('x',19)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULANDO CORRELAÇAO ENTRE AS SERIES  %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Vu - Vento Zonal
%  Vv - Vento Meridional
%  Cu - Corrente Zonal
%  Cv - Corrente Meridional
%  M - Mare
%  CRu - Corrente Residual Zonal
%  CRu - Corrente Residual Meridional
%  ElR - Elevacao Residual


disp('Correlacao Entre Vento Zonal - Corrente Zonal')
CC_Vu_Cu=corrcoef(estmet_u,corr_u);CC_Vu_Cu=abs(CC_Vu_Cu(2))*100

disp('Correlacao Entre Vento Merdional - Corrente Merdional')
CC_Vv_Cv=corrcoef(estmet_v,corr_v);CC_Vu_Cu=abs(CC_Vv_Cv(2))*100

disp('Correlacao Entre Vento Zonal - Nivel do Mar Residual')
CC_Vu_M=corrcoef(estmet_u,mare_amp);CC_Vu_M=abs(CC_Vu_M(2))*100

disp('Correlacao Entre Vento Merdional - Nivel do Mar Residual')
CC_Vv_M=corrcoef(estmet_v,mare_amp);CC_Vv_M=abs(CC_Vv_M(2))*100

disp('Correlacao Entre Corrente Zonal - Nivel do Mar Residual')
CC_Cu_M=corrcoef(corr_u,mare_amp);CC_Cu_M=abs(CC_Cu_M(2))*100

disp('Correlacao Entre Corrente Merdional - Nivel do Mar Residual Residual')
CC_Cv_M=corrcoef(corr_v,mare_amp);CC_Cv_M=abs(CC_Cv_M(2))*100


%  ------------------------

disp('Correlacao Entre Vento Zonal - Corrente Residual Zonal ')
CC_Vu_CRu=corrcoef(estmet_u,CRu);CC_Vu_CRu=abs(CC_Vu_CRu(2))*100

disp('Correlacao Entre Vento Merdional - Corrente Residual Merdional')
CC_Vv_CRv=corrcoef(estmet_v,CRv);CC_Vu_CRu=abs(CC_Vv_CRv(2))*100

disp('Correlacao Entre Vento Zonal - Nivel do Mar')
CC_Vu_ElR=corrcoef(estmet_u,ElR);CC_Vu_ElR=abs(CC_Vu_ElR(2))*100

disp('Correlacao Entre Vento Merdional - Nivel do Mar')
CC_Vv_ElR=corrcoef(estmet_v,ElR);CC_Vv_ElR=abs(CC_Vv_ElR(2))*100

disp('Correlacao Entre Corrente Zonal - Nivel do Mar')
CC_Cu_ElR=corrcoef(corr_u,ElR);CC_Cu_ElR=abs(CC_Cu_ElR(2))*100

disp('Correlacao Entre Corrente Merdional - Nivel do Mar')
CC_Cv_ElR=corrcoef(corr_v,ElR);CC_Cv_ElR=abs(CC_Cv_ElR(2))*100


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             AUTO-ESPECTROS		 %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=DATAC;

%  Frequencia
Fs=1/.5;
select = (1:(length(t)+1)/2)';
freq = (select - 1)*Fs/length(t);

%%%%%%%%%%%%%  CORRENTE ZONAL   %%%%%%%%%%%%%%

% faz a fft de y
fyCu=fft(corr_u); 

% modulo, elimina a fase (parte imag)
fyCu=fyCu.*conj(fyCu); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyCu=fliplr(fyCu((length(t)/2)+1:length(t)));

figure;  set(gcf,'color','w')
subplot(211)
plot(DATAC,corr_u)
title('Corrente Zonal [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyCu,'r');grid
title('Espectro da Corrente Zonal')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos



%%%%%%   CORRENTE MERIDIONAL  %%%%%%

% faz a fft de y
fyCv=fft(corr_v); 

% modulo, elimina a fase (parte imag)
fyCv=fyCv.*conj(fyCv); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyCv=fliplr(fyCv((length(t)/2)+1:length(t)));

figure;  set(gcf,'color','w')
subplot(211)
plot(DATAC,corr_v)
title('Corrente Meridional [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyCu,'r');grid
title('Espectro da Corrente Meridional')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos



%%%%%%%%  ELEVACAO DO NIVEL MAR  %%%%%%%%

% faz a fft de y
fyM=fft(mare_amp); 

% modulo, elimina a fase (parte imag)
fyM=fyM.*conj(fyM); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyM=fliplr(fyM((length(t)/2)+1:length(t)));

figure;  set(gcf,'color','w')
subplot(211)
plot(DATAC,mare_amp)
title('Elevacao do Nivel do Mar [m]')
xlabel('Data')
ylabel('Elevacao do Nivel do Mar [m]')
datetick('x',19)

subplot(212)
loglog(freq,fyM,'r');grid
title('Espectro da Elevacao do Nivel do Mar')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos


%%%%%%%%  VENTO ZONAL  %%%%%%%%

% faz a fft de y
fyVu=fft(estmet_u); 

% modulo, elimina a fase (parte imag)
fyVu=fyVu.*conj(fyVu); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyVu=fliplr(fyVu((length(t)/2)+1:length(t)));

figure;  set(gcf,'color','w')
subplot(211)
plot(DATAC,estmet_u)
title('Vento Zonal [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyVu,'r');grid
title('Espectro do Vento Zonal')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos


%%%%%%%%  VENTO MERIDIONAL  %%%%%%%%

% faz a fft de y
fyVv=fft(estmet_v); 

% modulo, elimina a fase (parte imag)
fyVv=fyVv.*conj(fyVv); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyVv=fliplr(fyVv((length(t)/2)+1:length(t)));

figure;    set(gcf,'color','w')
subplot(211)
plot(DATAC,estmet_v)
title('Vento Meridional [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyVv,'r');grid
title('Espectro do Vento Meridional')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos


%%%%%%%%%%%%%  CORRENTE RESIDUAL ZONAL   %%%%%%%%%%%%%%

% faz a fft de y
fyCRu=fft(CRu); 

% modulo, elimina a fase (parte imag)
fyCRu=fyCRu.*conj(fyCRu); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyCRu=fliplr(fyCRu((length(t)/2)+1:length(t)));

figure;   set(gcf,'color','w')
subplot(211)
plot(DATAC,CRu)
title('Corrente Residual Zonal [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyCu,'r');grid
title('Espectro da Corrente Residual Zonal')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos



%%%%%%   CORRENTE RESIDUAL MERIDIONAL  %%%%%%

% faz a fft de y
fyCRv=fft(CRv); 

% modulo, elimina a fase (parte imag)
fyCRv=fyCRv.*conj(fyCRv); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyCRv=fliplr(fyCRv((length(t)/2)+1:length(t)));

figure;   set(gcf,'color','w')
subplot(211)
plot(DATAC,CRv)
title('Corrente Residual Meridional [m.s^{-1}]')
xlabel('Data')
ylabel('Velocidade [m.s^{-1}]')
datetick('x',19)

subplot(212)
loglog(freq,fyCRv,'r');grid
title('Espectro da Corrente Residual Meridional')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos



%%%%%%%%  ELEVACAO DO NIVEL MAR RESIDUAL  %%%%%%%%

% faz a fft de y
fyElR=fft(ElR); 

% modulo, elimina a fase (parte imag)
fyElR=fyElR.*conj(fyElR); 

% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fyElR=fliplr(fyElR((length(t)/2)+1:length(t)));

figure;  set(gcf,'color','w')
subplot(211)
plot(DATAC,ElR)
title('Elevacao do Nivel do Mar  Residual [m]')
xlabel('Data')
ylabel('Elevacao do Nivel do Mar  Residual [m]')
datetick('x',19)

subplot(212)
loglog(freq,fyElR,'r');grid
title('Espectro da Elevacao do Nivel do Mar Residual')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

picos=ginput;
picos=picos(:,1);
picos=1./picos




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          ESPECTROS CRUZADOS		 %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  %  %  Vento Zonal - Corrente Zonal	%  %  %


% modulo, elimina a fase (parte imag)
fy_Vu_Cu=fyVu.*conj(fyCu);


% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fy_Vu_Cu=fliplr(fy_Vu_Cu((length(t)/2)+1:length(t))); 

% amplitude
xamp_Vu_Cu=sqrt(fy_Vu_Cu.*conj(fy_Vu_Cu));

%diagrama de fase
xphase_Vu_Cu=atan2(imag(fy_Vu_Cu),real(fy_Vu_Cu));

figure;  set(gcf,'color','w')

loglog(freq,xamp_Vu_Cu,'r')
hold on;  
axis('tight')
grid on
title('Espectro Cruzado do Vento e Corrente Zonais')
xlabel('Frequencia [horas]')
ylabel('Potencia Espectral')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  FAZER FASE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  %  %  Vento Merdional - Corrente Merdional	%  %  %




%  %  %  Vento Zonal - Nivel do Mar Residual	%  %  %




%  %  %  Vento Merdional - Nivel do Mar Residual	%  %  %




%  %  %  Corrente Zonal - Nivel do Mar Residual	   %  %  %




%  %  %  Corrente Merdional - Nivel do Mar Residual Residual	%  %  %




%  %  %  Vento Zonal - Corrente Residual Zonal	%  %  %




%  %  %  Vento Merdional - Corrente Residual Merdional	%  %  %




%  %  %  Vento Zonal - Nivel do Mar	%  %  %




%  %  %  Vento Merdional - Nivel do Mar	%  %  %




%  %  %  Corrente Zonal - Nivel do Mar	%  %  %




%  %  %  Corrente Merdional - Nivel do Mar	%  %  %










%  %%   Metodo 2 - PSD
%  
%  s=corr_u;
%  %s=MS;
%  % frequencia amostral
%  fs=1/(.5/24);
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  wind=ones(1,nfft)/nfft;
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,fs,wind,over,0.95);
%  
%  loglog(F,pxx);
%  
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cpd)^{-1} ]');
%  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cpd ] ');
%  set(hx,'FontSize',12,'FontWeight','demi')
%  grid on

