clear all;close all;clc
format short g

load corrF.mat
load corrS.mat
load mare.mat
load estmet.mat
load rot.mat

%  Corrente Residual
CRFu=corrF_u-pout_corrF_u';
CRFv=corrF_v-pout_corrF_v';

CRSu=corrS_u-pout_corrS_u';
CRSv=corrS_v-pout_corrS_v';

%  Elevacao Residual
ElR=mare_amp-pout_mare';


%  Vetor de DATA Completo
DATAC=datenum('23-Mar-2006 13:30'):datenum('00:30:00'):datenum('30-Apr-2006 23:31:00');

DATAM=datenum('17-Feb-2006 19:00'):datenum('00:30:00'):datenum('23-Mar-2006 10:30');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             AUTO-ESPECTROS		 %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=DATAC;

%  Frequencia


Fs=1/.5;
select = (1:(length(t)+1)/2)';
freq = (select - 1)*Fs/length(t);


%  Divisao da serie
nfft=floor(length(DATAC)/1);

%  Janela do filtro
jan=length(DATAC)/1;
%  jan=24;

%  tamanho do filtro
npf=24; % 12 horas
% tira o q esta em escesso na serie filtrada
inxs=(npf-1)/2; 

myfilter=hamming(npf);

myfilter=myfilter/sum(myfilter); 


%%%%%%   CORRENTE RESIDUAL SUPERFICIAL ZONAL  %%%%%%
%  
%  s=CRSu;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;				    % overlap=0
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);         % psd--> densidade de potencial espectral
%  						    % [Pxx, Pxxc, F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P)
%  						    % onde p é o intervalo de confiança
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRSu)
%  title('Corrente Residual Superficial Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual Superficial Zonal')
%  grid on
%  
%  picos=ginput;
%  picos=picos(:,1);
%  picos=1./picos
%  
%  %   print -depsc fyCRSu
%  %  !epstopdf fyCRSu.eps
%  
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT

%  % faz a fft de y
%  fyCRSu=fft(CRSu);
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRSu=fyCRSu.*conj(fyCRSu);
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRSu=fliplr(fyCRSu((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRSu)
%  title('Corrente Residual Superficial Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRSu,'r')
%  grid 
%  title('Espectro da Corrente Residual Superficial Zonal ')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRSu
%  !epstopdf fyCRSu.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos





%%%%%%   CORRENTE RESIDUAL SUPERFICIAL MERIDIONAL  %%%%%%
%  
%  s=CRSv;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRSv)
%  title('Corrente Residual Superficial Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  
%  subplot(212)
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
%  hy=ylabel('Densidade Espectro de Potencia [W.(kg.cph)^{-1}]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual Superficial Meridional')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono


%  % faz a fft de y
%  fyCRSv=fft(CRSv);
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRSv=fyCRSv.*conj(fyCRSv); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRSv=fliplr(fyCRSv((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRSv)
%  title('Corrente Residual Superficial Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRSv,'r')
%  grid
%  title('Espectro da Corrente Residual Superficial Meridional ')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRSv
%  !epstopdf fyCRSv.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos




%%%%%%   CORRENTE RESIDUAL DE FUNDO ZONAL  %%%%%%

%  s=CRFu;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRFu)
%  title('Corrente Residual de Fundo Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual de Fundo Zonal')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT  %%%%%

%  % faz a fft de y
%  fyCRFu=fft(CRFu);
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRFu=fyCRFu.*conj(fyCRFu); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRFu=fliplr(fyCRFu((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRFu)
%  title('Corrente Residual de Fundo Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRFu,'r');grid
%  title('Espectro da Corrente Residual de Fundo Zonal ')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRFu
%  !epstopdf fyCRFu.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos




%%%%%%   CORRENTE RESIDUAL DE FUNDO MERIDIONAL  %%%%%%
%  
%  s=CRFv;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRFv)
%  title('Corrente Residual de Fundo Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual de Fundo Meridional')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT  %%%%%

%  % faz a fft de y
%  fyCRFv=fft(CRFv);           
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRFv=fyCRFv.*conj(fyCRFv); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRFv=fliplr(fyCRFv((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRFv)
%  title('Corrente Residual de Fundo Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRFv,'r');grid
%  title('Espectro da Corrente Residual de Fundo Meridional ')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRFv
%  !epstopdf fyCRFv.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos



%%%%%%   CORRENTE RESIDUAL SUPERFICIAL ROTACIONADA ALONG CHANNEL  %%%%%%

s=CRS_rot_ac;
y=conv(s,myfilter); 
% remove os excessos das pontas 
y=y(inxs+1:length(y)-inxs); 



wind=ones(1,jan); %/nfft;
%  wind=hamming(npf);
over=0; %over=nfft/2;
[pxx,pxxc,F]=psd(y,nfft,Fs,wind,over,0.68);

figure
%  subplot(211)
%  plot(DATAC,CRS_rot_ac);hold on,plot(DATAC,y,'r')
%  title('Corrente Residual Superficial Rotacionada Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on

%  subplot(212)
loglog(F,pxx);
% plota intervalo de confianca
F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];

hold on
hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
hp=loglog(F,pxx,'k');
set(hp,'LineWidth',2.0)
set(hpa,'EdgeColor',[1 1 1])
hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  set(hy,'FontSize',12,'FontWeight','demi')
hx=xlabel(' Frequencia [ cph ] ');
%  set(hx,'FontSize',12,'FontWeight','demi')
title('Espectro da Corrente Residual Superficial Rotacionada Along Channel')
grid on

print -depsc Cs_ac
!epstopdf Cs_ac.eps
!mv Cs_ac.pdf figuras/.


%  picos=ginput;
%  picos=picos(:,1);
%  picos=1./picos

clear s F pxx pxxc poligono y


%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyCRS_rot_ac=fft(CRS_rot_ac);           
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRS_rot_ac=fyCRS_rot_ac.*conj(fyCRS_rot_ac); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRS_rot_ac=fliplr(fyCRS_rot_ac((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRS_rot_ac)
%  title('Corrente Residual Superficial Rotacionada Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRS_rot_ac,'r');grid
%  title('Espectro da Corrente Residual Superficial Rotacionada Along Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRS_rot_ac
%  !epstopdf fyCRS_rot_ac.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos




%%%%%%   CORRENTE RESIDUAL SUPERFICIAL ROTACIONADA CROSS CHANNEL  %%%%%%
%  
%  s=CRS_rot_cc;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRS_rot_cc)
%  title('Corrente Residual Superficial Rotacionada Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual Superficial Rotacionada Cross Channel')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyCRS_rot_cc=fft(CRS_rot_cc);           
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRS_rot_cc=fyCRS_rot_cc.*conj(fyCRS_rot_cc); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRS_rot_cc=fliplr(fyCRS_rot_cc((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRS_rot_cc)
%  title('Corrente Residual Superficial Rotacionada Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRS_rot_cc,'r');grid
%  title('Espectro da Corrente Residual Superficial Rotacionada Cross Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRS_rot_cc
%  !epstopdf fyCRS_rot_cc.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  



%%%%%%   CORRENTE RESIDUAL DE FUNDO ROTACIONADA ALONG CHANNEL  %%%%%%

s=CRF_rot_ac;
y=conv(s,myfilter); 
% remove os excessos das pontas 
y=y(inxs+1:length(y)-inxs); 


wind=ones(1,jan); %/nfft;
%  wind=hamming(npf);
over=0; %over=nfft/2;
[pxx,pxxc,F]=psd(y,nfft,Fs,wind,over,0.68);

figure
%  subplot(211)
%  plot(DATAC,CRF_rot_ac);hold on,plot(DATAC,y,'r')
%  title('Corrente Residual de Fundo Rotacionada Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on

%  subplot(212)
loglog(F,pxx);
% plota intervalo de confianca
F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];

hold on
hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
hp=loglog(F,pxx,'k');
set(hp,'LineWidth',2.0)
set(hpa,'EdgeColor',[1 1 1])
hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  set(hy,'FontSize',12,'FontWeight','demi')
hx=xlabel(' Frequencia [ cph ] ');
%  set(hx,'FontSize',12,'FontWeight','demi')
title('Espectro da Corrente de Fundo Superficial Rotacionada Along Channel')
grid on

print -depsc Cf_ac
!epstopdf Cf_ac.eps
!mv Cf_ac.pdf figuras/.

%  picos=ginput;
%  picos=picos(:,1);
%  picos=1./picos

clear s F pxx pxxc poligono y 

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyCRF_rot_ac=fft(CRF_rot_ac);
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRF_rot_ac=fyCRF_rot_ac.*conj(fyCRF_rot_ac); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRF_rot_ac=fliplr(fyCRF_rot_ac((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRF_rot_ac)
%  title('Corrente Residual de Fundo Rotacionada Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRF_rot_ac,'k');grid
%  title('Espectro da Corrente Residual de Fundo Rotacionada Along Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRF_rot_ac
%  !epstopdf fyCRF_rot_ac.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos





%%%%%%   CORRENTE RESIDUAL DE FUNDO ROTACIONADA CROSS CHANNEL  %%%%%%
%  
%  s=CRF_rot_cc;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRF_rot_cc)
%  title('Corrente Residual de Fundo Rotacionada Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Corrente Residual de Fundo Rotacionada Cross Channel')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyCRF_rot_cc=fft(CRF_rot_cc);           
%  
%  % modulo, elimina a fase (parte imag)
%  fyCRF_rot_cc=fyCRF_rot_cc.*conj(fyCRF_rot_cc); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyCRF_rot_cc=fliplr(fyCRF_rot_cc((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,CRF_rot_cc)
%  title('Corrente Residual de Fundo Rotacionada Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyCRF_rot_cc,'k');grid
%  title('Espectro da Corrente Residual de Fundo Rotacionada Cross Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyCRF_rot_cc
%  !epstopdf fyCRF_rot_cc.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos


%%%%%%%%  VENTO ZONAL  %%%%%%%%
%  
%  s=estmet_u;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,estmet_u)
%  title('Vento Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro do Vento Zonal')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT

%  % faz a fft de y
%  fyVu=fft(estmet_u); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyVu=fyVu.*conj(fyVu); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyVu=fliplr(fyVu((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,estmet_u)
%  title('Vento Zonal [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyVu,'m');grid
%  title('Espectro do Vento Zonal')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyVu
%  !epstopdf fyVu.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos



%%%%%%%%  VENTO MERIDIONAL  %%%%%%%%
%  s=estmet_v;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,estmet_v)
%  title('Vento Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro do Vento Meridional')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyVv=fft(estmet_v); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyVv=fyVv.*conj(fyVv); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyVv=fliplr(fyVv((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,estmet_v)
%  title('Vento Meridional [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyVv,'m');grid
%  title('Espectro do Vento Meridional')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyVv
%  !epstopdf fyVv.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos



%%%%%%%%  VENTO ROTACIONADO ALONG CHANNEL  %%%%%%%%
s=V_rot_ac;
y=conv(s,myfilter); 
% remove os excessos das pontas 
y=y(inxs+1:length(y)-inxs); 



wind=ones(1,jan); %/nfft;
%  wind=hamming(npf);
over=0; %over=nfft/2;
[pxx,pxxc,F]=psd(y,nfft,Fs,wind,over,0.68);

figure
%  subplot(211)
%  plot(DATAC,V_rot_ac);hold on,plot(DATAC,y,'r')
%  title('Vento Rotacionado Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on

%  subplot(212)
loglog(F,pxx);
% plota intervalo de confianca
F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];

hold on
hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
hp=loglog(F,pxx,'k');
set(hp,'LineWidth',2.0)
set(hpa,'EdgeColor',[1 1 1])
hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  set(hy,'FontSize',12,'FontWeight','demi')
hx=xlabel(' Frequencia [ cph ] ');
%  set(hx,'FontSize',12,'FontWeight','demi')
title('Espectro do Vento Rotacionado Along Channel')
grid on

%  picos=ginput;
%  picos=picos(:,1);
%  picos=1./picos


print -depsc V_ac
!epstopdf V_ac.eps
!mv V_ac.pdf figuras/.

clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyV_rot_ac=fft(V_rot_ac); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyV_rot_ac=fyV_rot_ac.*conj(fyV_rot_ac); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyV_rot_ac=fliplr(fyV_rot_ac((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,V_rot_ac)
%  title('Vento Rotacionado Along Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyV_rot_ac,'m');grid
%  title('Espectro do Vento Rotacionado Along Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyV_rot_ac
%  !epstopdf fyV_rot_ac.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos


%%%%%%%%  VENTO ROTACIONADO CROSS CHANNEL  %%%%%%%%
%  
%  s=V_rot_cc;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAC,V_rot_cc)
%  title('Vento Rotacionado Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro do Vento Rotacionado Cross Channel')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyV_rot_cc=fft(V_rot_cc); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyV_rot_cc=fyV_rot_cc.*conj(fyV_rot_cc); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyV_rot_cc=fliplr(fyV_rot_cc((length(t)/2)+1:length(t)));
%  
%  figure
%  subplot(211)
%  plot(DATAC,V_rot_cc)
%  title('Vento Rotacionado Cross Channel [m.s^{-1}]')
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freq(1:end-1),fyV_rot_cc,'m');grid
%  title('Espectro do Vento Rotacionado Cross Channel')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyV_rot_cc
%  !epstopdf fyV_rot_cc.eps
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos




%  %%%%%%%%  ELEVACAO DO NIVEL MAR  %%%%%%%%

%  tM=DATAM;
%  
%  %  Frequencia
%  selectM = (1:(length(tM)+1)/2)';
%  freqM = (selectM - 1)*Fs/length(tM);
%  
%  s=mare_amp;
%  
%  nfft=floor(length(s)/1);
%  % wind=window(@hamming,floor(nfft/2));
%  %wind=window(@hamming,nfft);
%  jan=length(s)/1;
%  
%  wind=ones(1,jan); %/nfft;
%  %  wind=hamming(npf);
%  over=0; %over=nfft/2;
%  [pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);
%  
%  figure
%  subplot(211)
%  plot(DATAM,mare_amp)
%  title('Elevação do Nível do Mar [m.s^{-1}]')         % ATENÇÃO: mudar a unidade
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
%  loglog(F,pxx);
%  % plota intervalo de confianca
%  F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
%  poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];
%  
%  hold on
%  hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
%  hp=loglog(F,pxx,'k');
%  set(hp,'LineWidth',2.0)
%  set(hpa,'EdgeColor',[1 1 1])
%  hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  %  set(hy,'FontSize',12,'FontWeight','demi')
%  hx=xlabel(' Frequencia [ cph ] ');
%  %  set(hx,'FontSize',12,'FontWeight','demi')
%  title('Espectro da Elevação do Nível do Mar')
%  grid on
%  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos
%  
%  clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%
  
%  % faz a fft de y
%  fyM=fft(mare_amp); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyM=fyM.*conj(fyM); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyM=fliplr(fyM((length(tM)/2)+1:length(tM)));
%  
%  figure
%  subplot(211)
%  plot(DATAM,mare_amp)
%  title('Elevacao do Nivel do Mar [m]')
%  xlabel('Data')
%  ylabel('Elevacao do Nivel do Mar [m]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freqM,fyM,'b');grid
%  title('Espectro da Elevacao do Nivel do Mar')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyM
%  !epstopdf fyM.eps
%  %  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos


%  %%%%%%%%  ELEVACAO DO NIVEL DO MAR RESIDUAL  %%%%%%%%

s=ElR;
%  Divisao da serie
nfft=floor(length(s)/1);

%  Janela do filtro
jan=length(s)/1;
%  jan=24;

wind=ones(1,jan); %/nfft;
%  wind=hamming(npf);
over=0; %over=nfft/2;
[pxx,pxxc,F]=psd(s,nfft,Fs,wind,over,0.68);

figure
%  subplot(211)
%  plot(DATAM,ElR)
%  title('Elevacao do Nivel do Mar Residual [m.s^{-1}]')         % ATENÇÃO: mudar a unidade
%  xlabel('Data')
%  ylabel('Velocidade [m.s^{-1}]')
%  datetick('x',19);axis tight;grid on
%  
%  subplot(212)
loglog(F,pxx);
% plota intervalo de confianca
F=F(2:end); pxxc=pxxc(2:end,:);pxx=pxx(2:end);
poligono=[F pxxc(:,1);flipud(F) flipud(pxxc(:,2));F(1)    pxxc(1,1)];

hold on
hpa= patch(poligono(:,1), poligono(:,2),[0 0.8 0.8]);
hp=loglog(F,pxx,'k');
set(hp,'LineWidth',2.0)
set(hpa,'EdgeColor',[1 1 1])
hy=ylabel('Densidade Espectro de Potencia [ W (kg.cph)^{-1} ]');
%  set(hy,'FontSize',12,'FontWeight','demi')
hx=xlabel(' Frequencia [ cph ] ');
%  set(hx,'FontSize',12,'FontWeight','demi')
title('Espectro da Elevacao do Nivel do Mar Residual')
grid on

print -depsc ElR
!epstopdf ElR.eps
!mv ElR.pdf figuras/.

%  picos=ginput;
%  picos=picos(:,1);
%  picos=1./picos

clear s F pxx pxxc poligono

%%%%  METODO FFT %%%%%%

%  % faz a fft de y
%  fyElR=fft(ElR); 
%  
%  % modulo, elimina a fase (parte imag)
%  fyElR=fyElR.*conj(fyElR); 
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fyElR=fliplr(fyElR((length(tM)/2)+1:length(tM)));
%  
%  figure
%  subplot(211)
%  plot(DATAM,ElR)
%  title('Elevacao do Nivel do Mar  Residual [m]')
%  xlabel('Data')
%  ylabel('Elevacao do Nivel do Mar  Residual [m]')
%  datetick('x',19)
%  
%  subplot(212)
%  loglog(freqM,fyElR,'b');grid
%  title('Espectro da Elevacao do Nivel do Mar Residual')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fyElR
%  !epstopdf fyElR.eps
%  %  
%  %  picos=ginput;
%  %  picos=picos(:,1);
%  %  picos=1./picos

