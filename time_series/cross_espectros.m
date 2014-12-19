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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          ESPECTROS CRUZADOS		 %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=DATAC;

%  Frequencia
Fs=1/.5;
select = (1:(length(t)+1)/2)';
freq = (select - 1)*Fs/length(t);

%  tamanho do filtro
npf=24; % 12 horas
% tira o q esta em escesso na serie filtrada
inxs=(npf-1)/2; 

myfilter=hamming(npf);

myfilter=myfilter/sum(myfilter); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  %  %  Vento Rotacionado Along Channel - Corrente  Residual Superficial Rotacionada Along Channel	%  %  %
s1=V_rot_ac;
s2=CRS_rot_ac;
% faz a fft de V_rot_ac e CRS_rot_ac
y1=conv(s1,myfilter); 
% remove os excessos das pontas 
y1=y1(inxs+1:length(y1)-inxs); 

y2=conv(s2,myfilter); 
% remove os excessos das pontas 
y2=y2(inxs+1:length(y2)-inxs); 


fyV_rot_ac=fft(y1); 
fyCRS_rot_ac=fft(y2); 

% modulo, elimina a fase (parte imag)
fy_V_rot_ac_CRS_rot_ac=fyV_rot_ac'.*conj(fyCRS_rot_ac);


% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fy_V_rot_ac_CRS_rot_ac=fliplr(fy_V_rot_ac_CRS_rot_ac((length(t)/2)+1:length(t))); 


% amplitude
xamp_V_rot_ac_CRS_rot_ac=sqrt(fy_V_rot_ac_CRS_rot_ac.*conj(fy_V_rot_ac_CRS_rot_ac));

%diagrama de fase
xphase_V_rot_ac_CRS_rot_ac=atan2(imag(fy_V_rot_ac_CRS_rot_ac),real(fy_V_rot_ac_CRS_rot_ac));

figure
loglog(freq(1:end-1),xamp_V_rot_ac_CRS_rot_ac,'b')
axis('tight')
grid on
title('Espectro Cruzado do Vento Along Channel e Corrente Residual Superficial Along Channel ')
xlabel('Frequencia [cph]')
ylabel('Potencia Espectral')

print -depsc fy_V_rot_ac_CRS_rot_ac
!epstopdf fy_V_rot_ac_CRS_rot_ac.eps
%  !mv fy_V_rot_ac_CRS_rot_ac.pdf figuras/.

%  loglog(freq(38),xamp_V_rot_ac_CRS_rot_ac(38),'r*')

figure
semilogx(freq(1:end-1),(180/pi)*xphase_V_rot_ac_CRS_rot_ac,'r');hold on
semilogx(freq(38),(180/pi)*xphase_V_rot_ac_CRS_rot_ac(38),'bo');hold on 
plot([freq(38) freq(38)],[-180 180],'b'); grid on
axis('tight')
title('Fase do Espectro Cruzado entre o Vento Along Channel e aCorrente Residual Superficial Along Channel ')
xlabel('Frequencia [cph]')
ylabel('Angulo')

print -depsc fase_fy_V_rot_ac_FCRS_rot_ac
!epstopdf fase_fy_V_rot_ac_FCRS_rot_ac.eps


stop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  %  Vento Rotacionado Along Channel - Corrente  Residual Fundo Rotacionada Along Channel	%  %  %
s1=V_rot_ac;
s2=CRF_rot_ac;
% faz a fft de V_rot_ac e CRS_rot_ac
y1=conv(s1,myfilter); 
% remove os excessos das pontas 
y1=y1(inxs+1:length(y1)-inxs); 

y2=conv(s2,myfilter); 
% remove os excessos das pontas 
y2=y2(inxs+1:length(y2)-inxs); 



% faz a fft de V_rot_ac e CRS_rot_ac
fyV_rot_ac=fft(y1); 
fyCRF_rot_ac=fft(y2); 

% modulo, elimina a fase (parte imag)
fy_V_rot_ac_CRF_rot_ac=fyV_rot_ac'.*conj(fyCRF_rot_ac);


% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fy_V_rot_ac_CRF_rot_ac=fliplr(fy_V_rot_ac_CRF_rot_ac((length(t)/2)+1:length(t))); 


% amplitude
xamp_V_rot_ac_CRF_rot_ac=sqrt(fy_V_rot_ac_CRF_rot_ac.*conj(fy_V_rot_ac_CRF_rot_ac));

%diagrama de fase
xphase_V_rot_ac_CRF_rot_ac=atan2(imag(fy_V_rot_ac_CRF_rot_ac),real(fy_V_rot_ac_CRF_rot_ac));

figure
loglog(freq(1:end-1),xamp_V_rot_ac_CRF_rot_ac,'b')
axis('tight')
grid on
title('Espectro Cruzado do Vento Along Channel e Corrente Residual Fundo Along Channel ')
xlabel('Frequencia [cph]')
ylabel('Potencia Espectral')

print -depsc fy_V_rot_ac_CRF_rot_ac
!epstopdf fy_V_rot_ac_CRF_rot_ac.eps
%  !mv fy_V_rot_ac_CRF_rot_ac.pdf figuras/.


%  loglog(freq(38),xamp_V_rot_ac_CRF_rot_ac(38),'r*'


figure
semilogx(freq(1:end-1),(180/pi)*xphase_V_rot_ac_CRF_rot_ac,'r');hold on
semilogx(freq(38),(180/pi)*xphase_V_rot_ac_CRF_rot_ac(38),'bo');hold on 
plot([freq(38) freq(38)],[-180 180],'b'); grid on
axis('tight')
grid on
title('Fase do Espectro Cruzado entre o Vento Along Channel e aCorrente Residual Fundo Along Channel ')
xlabel('Frequencia [cph]')
ylabel('Angulo')

print -depsc fase_fy_V_rot_ac_CRF_rot_ac
!epstopdf fase_fy_V_rot_ac_CRF_rot_ac.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %  %  %  %  Corrente Residual Superficial Along Channel Zonal  - Corrente Residual de Fundo Along Channel   %  %  %
%  

s1=CRS_rot_ac;
s2=CRF_rot_ac;
% faz a fft de V_rot_ac e CRS_rot_ac
y1=conv(s1,myfilter); 
% remove os excessos das pontas 
y1=y1(inxs+1:length(y1)-inxs); 

y2=conv(s2,myfilter); 
% remove os excessos das pontas 
y2=y2(inxs+1:length(y2)-inxs); 


% faz a fft de CRSu e CRFu
fyCRSu=fft(y1); 
fyCRFu=fft(y2);

% modulo, elimina a fase (parte imag)
fy_CRSu_CRFu=fyCRSu.*conj(fyCRFu);


% o comando fliplr altera a ordem da colunas 
% sendo simetrica, pegamos 1/2 de um lado
fy_CRSu_CRFu=fliplr(fy_CRSu_CRFu((length(t)/2)+1:length(t))); 


% amplitude
xamp_CRSu_CRFu=sqrt(fy_CRSu_CRFu.*conj(fy_CRSu_CRFu));

%diagrama de fase
xphase_CRu_CRFu=atan2(imag(fy_CRSu_CRFu),real(fy_CRSu_CRFu));

figure
loglog(freq(1:end-1),xamp_CRSu_CRFu,'b')
axis('tight')
grid on
title('Espectro Cruzado da Corrente Residual Superficial Along Channel e Corrente Residual de Fundo Along Channel')
xlabel('Frequencia [cph]')
ylabel('Potencia Espectral')
 
print -depsc fy_CRS_rot_CRF_rot
!epstopdf fy_CRS_rot_CRF_rot.eps
%  !mv fy_CRS_rot_CRF_rot.pdf figuras/.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%  FAZER FASE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%  
%  %  %  %  Vento Rotacionado Cross Channel - Corrente  Residual Superficial Rotacionada Cross Channel	%  %  %
%  
%  % faz a fft de V_rot_ac e CRS_rot_ac
%  fyV_rot_cc=fft(V_rot_cc); 
%  fyCRS_rot_cc=fft(CRS_rot_cc); 
%  
%  % modulo, elimina a fase (parte imag)
%  fy_V_rot_cc_CRS_rot_cc=fyV_rot_cc'.*conj(fyCRS_rot_cc);
%  
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fy_V_rot_cc_CRS_rot_cc=fliplr(fy_V_rot_cc_CRS_rot_cc((length(t)/2)+1:length(t))); 
%  
%  
%  % amplitude
%  xamp_V_rot_cc_CRS_rot_cc=sqrt(fy_V_rot_cc_CRS_rot_cc.*conj(fy_V_rot_cc_CRS_rot_cc));
%  
%  %diagrama de fase
%  xphase_V_rot_cc_CRS_rot_cc=atan2(imag(fy_V_rot_cc_CRS_rot_cc),real(fy_V_rot_cc_CRS_rot_cc));
%  
%  figure
%  loglog(freq(1:end-1),xamp_V_rot_cc_CRS_rot_cc,'b')
%  axis('tight')
%  grid on
%  title('Espectro Cruzado do Vento Rotacionado Cross Channel e Corrente Residual Superficial Rotacionada Cross Channel ')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%  
%  print -depsc fy_V_rot_cc_CRS_rot_cc
%  !epstopdf fy_V_rot_cc_CRS_rot_cc.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%  FAZER FASE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



%  
%  %  %  %  %  Corrente Residual Superficial Zonal  - Corrente Residual de Fundo Zonal   %  %  %
%  
%  % faz a fft de CRSu e CRFu
%  fyCRSu=fft(CRSu); 
%  fyCRFu=fft(CRFu);
%  
%  % modulo, elimina a fase (parte imag)
%  fy_CRSu_CRFu=fyCRSu.*conj(fyCRFu);
%  
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fy_CRSu_CRFu=fliplr(fy_CRSu_CRFu((length(t)/2)+1:length(t))); 
%  
%  
%  % amplitude
%  xamp_CRSu_CRFu=sqrt(fy_CRSu_CRFu.*conj(fy_CRSu_CRFu));
%  
%  %diagrama de fase
%  xphase_CRu_CRFu=atan2(imag(fy_CRSu_CRFu),real(fy_CRSu_CRFu));
%  
%  figure
%  loglog(freq(1:end-1),xamp_CRSu_CRFu,'g')
%  axis('tight')
%  grid on
%  title('Espectro Cruzado da Corrente Residual Superficial Zonal e Corrente Residual de Fundo Zonal')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%   
%  print -depsc fy_CRSu_CRFu
%  !epstopdf fy_CRSu_CRFu.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%  FAZER FASE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%  
%  %  %  %  %  Corrente Residual Superficial Meridional  - Corrente Residual de Fundo Meridional   %  %  %
%  
%  % faz a fft de CRSv e CRFv
%  fyCRSv=fft(CRSv); 
%  fyCRFv=fft(CRFv);
%  
%  % modulo, elimina a fase (parte imag)
%  fy_CRSv_CRFv=fyCRSv.*conj(fyCRFv);
%  
%  
%  % o comando fliplr altera a ordem da colunas 
%  % sendo simetrica, pegamos 1/2 de um lado
%  fy_CRSv_CRFv=fliplr(fy_CRSv_CRFv((length(t)/2)+1:length(t))); 
%  
%  
%  % amplitude
%  xamp_CRSv_CRFv=sqrt(fy_CRSv_CRFv.*conj(fy_CRSv_CRFv));
%  
%  %diagrama de fase
%  xphase_CRv_CRFv=atan2(imag(fy_CRSv_CRFv),real(fy_CRSv_CRFv));
%  
%  figure
%  loglog(freq(1:end-1),xamp_CRSv_CRFv,'g')
%  axis('tight')
%  grid on
%  title('Espectro Cruzado da Corrente Residual Superficial Meridional e Corrente Residual de Fundo Meridional')
%  xlabel('Frequencia [horas]')
%  ylabel('Potencia Espectral')
%   
%  print -depsc fy_CRSv_CRFv
%  !epstopdf fy_CRSv_CRFv.eps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%  FAZER FASE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



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

