%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELBORACAO DE FIGURA EVIDENCIANDO 
% AS ETAPAS DE TRATAMENTO BASICO
%        DE DADOS DE CTD
%
% Rafael Soutelino - Mestrado IOUSP
%           julho/2007
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

xx = 45;
XX = num2str(xx);

if xx >= 1 & xx <= 9
    fid=fopen(['/home/rafaelgs/mestrado/dados/leste2/ctd/brutos/OE00',XX,'.cnv'],'r');
    for i = 1:49; fgetl(fid); end
    ctd = fscanf(fid,'%f',[6 inf]);
    ctd = ctd';
  elseif  xx >= 10 & xx <= 99
    fid=fopen(['/home/rafaelgs/mestrado/dados/leste2/ctd/brutos/OE0',XX,'.cnv'],'r');
    for i = 1:49; fgetl(fid); end
    ctd = fscanf(fid,'%f',[6 inf]);
    ctd = ctd';
  else 
    fid=fopen(['/home/rafaelgs/mestrado/dados/leste2/ctd/brutos/OE',XX,'.cnv'],'r');
    for i = 1:49; fgetl(fid); end
    ctd = fscanf(fid,'%f',[6 inf]);
    ctd = ctd';
  end

  p=ctd(:,1);
  T=ctd(:,2);
  C=ctd(:,4);

   %% separando os dados de descida
  d = find(p == max(p));  
  p = p(1:d);
  T = T(1:d);
  C = C(1:d);
  C = C*1000/100; % convertendo para mS/cm

  %%% CONVERSAO DA CONDUTIVIDADE EM SALINIDADE %%%
  
  R=1/sw_c3515*C;
  S=sw_salt(R,T,p);

  figure(1)
  set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 13 8],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 13 8],...
        'ShareColors','off',...
        'Clipping','on');
  
  subplot(131)
  set(gca,'fontsize',14)
  plot(T,-p,'.');hold on
  axis([23 29 -120 -20])
  ylabel('Profundidade [m]','fontweight','bold')
  title('Bruto','fontweight','bold')
  pbaspect([0.5 0.8 0.5])

 %%% PREPARACAO PARA MEDIA EM CAIXA %%%
  
  maxp=max(p); minp=min(p);
  % motando caixas de 1 dbar
  pmin=ceil(minp)+.5;
  pmax=floor(maxp)-.5;

  f=find(p>=pmin & p<=pmax); 
  p=p(f);
  T=T(f);
  S=S(f);

  CC=C(f);
  l1=length(p);

  clear ctd

  binp=pmin:pmax;
  lp=length(binp);

  %%% MEDIA EM CAIXA %%%
  
  for kk=1:lp-1
    f=find(p>=binp(kk) & p<=binp(kk+1));
    Tbin(kk)=mean(T(f));
    CCbin(kk)=mean(CC(f));
    Sbin(kk)=mean(S(f));
    pbin(kk)=.5*(binp(kk)+binp(kk+1));
  end

  subplot(132)
  set(gca,'fontsize',14)
  plot(Tbin,-pbin,'.');hold on
  axis([23 29 -120 -20])
  pbaspect([0.5 0.8 0.5])
  title('Binado','fontweight','bold')
  xlabel('Temperatura [^\circ C]','fontweight','bold')

  %%% JANELA MOVEL %%%

    % basicamente o procedimento de convolucao
    Ts=weim(31,'hann',Tbin);
    CCs=weim(31,'hann',CCbin);
    Ss=weim(31,'hann',Sbin);

  subplot(133)
  set(gca,'fontsize',14)
  plot(Tbin,-pbin,'.');hold on
  plot(Ts,-pbin,'r','linewidth',2)
  axis([23 29 -120 -20])
  pbaspect([0.5 0.8 0.5])
  title('Filtrado','fontweight','bold')

print -depsc ../figuras/fig_trat_basico.eps
!epstopdf ../figuras/fig_trat_basico.eps






