%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                               
%        TRATAMENTO BASICO DE DADOS HIDROGRAFICOS
%                    OCEANO LESTE II
%                 out / 2006 - Mestrado
%                   Rafael Soutelino                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  - separa os dados de descida
%  - faz media em caixa
%  - aplica janela movel


clear;close all;clc;Tm=[];Sm=[];pm=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOPING PRINCIPAL DE PROCESSAMENTO DE CADA ESTACAO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx=input('Estacoes a serem tratadas [1:16]: ');
lx=length(xx);

%%% OBS: a rotina toda é um looping a partir daqui

for k=1:lx;
  XX=num2str(xx(k));
  
  if xx(k) >= 1 & xx(k) <= 9
    fid=fopen(['../../../dados/leste2/ctd/brutos/OE00',XX,'.cnv'],'r');
    for i = 1:49; fgetl(fid); end
    ctd = fscanf(fid,'%f',[6 inf]);
    ctd = ctd';
  elseif  xx(k) >= 10 & xx(k) <= 99
    fid=fopen(['../../../dados/leste2/ctd/brutos/OE0',XX,'.cnv'],'r');
    for i = 1:49; fgetl(fid); end
    ctd = fscanf(fid,'%f',[6 inf]);
    ctd = ctd';
  else 
    fid=fopen(['../../../dados/leste2/ctd/brutos/OE',XX,'.cnv'],'r');
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

  %%% JANELA MOVEL %%%

    % basicamente o procedimento de convolucao
    Ts=weim(31,'hann',Tbin);
    CCs=weim(31,'hann',CCbin);
    Ss=weim(31,'hann',Sbin);

 
  Tbru=T;Sbru=S;pbru=p; % renomeando os vetores para os perfis brutos
  T=Ts';S=Ss';p=pbin';C=CCs'; % é necessário transpor para salvar depois


  %%% CALCULO DE Sigma-T %%%
  sigT=sw_dens(S,T,p)-1000;

  %%% CALCULO DE Sigma-Theta %%%
  sigTH=sw_pden(S,T,p,0)-1000;

  
  %%% PLOTAGEM DE PERFIS %%%
  
    t1=['Estacao ',XX,' - Temperatura com \Deltaz = 1 m'];
    t2=['Estacao ',XX,' - Salinidade com \Deltaz = 1 m'];
    t3=['Estacao ',XX,' - Densidades \sigma_T e \sigma_\theta com \Deltaz = 1 m'];
  
    figure(2)
    plot(Tbru,-pbru,'b--','linewidth',.5)
    hold on
    plot(Ts,-pbin,'r','linewidth',1);
    hold off
    title(t1,'fontsize',10,'fontweight','bold')
    xlabel('Temperatura [ \circ C]','fontsize',10)
    ylabel('Profundidade (m)','fontsize',10)
    set(gca,'PlotBoxAspectRatio',[1 2 1],'TickLength',[0 0])
    drawnow
%    print(2,'-depsc',['figuras/perfil_temp_',XX]);
%    eval(['!epstopdf figuras/perfil_temp_',XX,'.eps'])
%    print(2,'-dpng',['../figuras/perfil_temp_',XX]);
%    
    figure(3)
    plot(Sbru,-pbru,'b--','linewidth',.5)
    hold on
    plot(Ss,-pbin,'r','linewidth',1);
    hold off
    title(t2,'fontsize',10,'fontweight','bold')
    xlabel('Salinidade','fontsize',10)
    ylabel('Profundidade (m)','fontsize',10)
    set(gca,'PlotBoxAspectRatio',[1 2 1],'TickLength',[0 0])
    drawnow
%    print(3,'-depsc',['figuras/perfil_sal_',XX]);
%    eval(['!epstopdf figuras/perfil_sal_',XX,'.eps'])
%    print(3,'-dpng',['../figuras/perfil_sal_',XX]);

%% nao fez diferença para coluna tao rasa
    figure(4)
    plot(sigT,-pbin,'g','linewidth',1);
    hold on
    plot(sigTH,-pbin,'m','linewidth',1);
    leg=legend('\sigma_T','\sigma_\theta',1);
    set(leg,'fontsize',12);
    hold off
    title(t3,'fontsize',10,'fontweight','bold')
    xlabel('Densidade [kg m^{-3}]','fontsize',10)
    ylabel('Profundidade (m)','fontsize',10)
    set(gca,'PlotBoxAspectRatio',[1 2 1],'TickLength',[0 0])
    drawnow
%    print(4,'-depsc',['figuras/perfil_dens_',XX]);
%    eval(['!epstopdf figuras/perfil_dens_',XX,'.eps']) 
%    print(4,'-dpng',['../figuras/perfil_dens_',XX]);

  
  %%% SALVANDO OS DADOS TRATADOS %%%   % comentado, pois os dados ja foram salvos


  format long g
  data=[p'; T'; S']; 

  eval(['fid=fopen(''../../../dados/leste2/ctd/filtrados/lesteII_ctd',XX,'.dat'',''w'');']);

  fprintf(fid,'%12.4f %12.4f %12.4f\n',data);

  clear p T S C CC CCbin CCs binp pbin Tbin Sbin Ss Ts pmax pmin maxp kk
  clear minpl1 lp lx t1 t2 t3 Tgarf Sgarf Siggarf garf Garf garftest
  clear R RR Sbru Tbru pbru Sg TTH Tg XX data minp q sigT sigTH svel
  
end
