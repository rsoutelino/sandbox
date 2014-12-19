%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                               
%        TRATAMENTO BASICO DE DADOS HIDROGRAFICOS
%                    SAO SEBASTIAO
%             maio / 2006 - Observacional                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  - separa os dados de descida
%  - faz media em caixa
%  - aplica janela movel
%  - salva os dados com o nome da radial que eles pertencem

clear;close all;clc;Tm=[];Sm=[];pm=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOPING PRINCIPAL DE PROCESSAMENTO DE CADA ESTACAO %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  xx=[1:28]; % para usar esse, trocar 70 por 72 e 10 por 11 no else
xx=[29:51 53:63 65 67 68 70:80 82:113]; % para usar esse, trocar 72 por 70 e 11 por 10 no else
lx=length(xx);

%%% OBS: a rotina toda é um looping a partir daqui

for k=1:lx;
  XX=num2str(xx(k));
  if xx(k)==1 
     fid=fopen(['ctd/brutos/doei0',XX,'.cnv'],'r');
     for i = 1:70; fgetl(fid); end
     ctd = fscanf(fid,'%f',[10 inf]);
     ctd = ctd';
  elseif xx(k)==19 | xx(k)==20 | xx(k)==21 | xx(k)==22 | xx(k)==23 | xx(k)==24 | xx(k)==25 | xx(k)==26 | xx(k)==27 | xx(k)==28 
     eval(['load ../../../dados/leste1/brutos/doei0',XX,'.cnv']);
     eval(['ctd = doei0',XX,';']);
     eval(['clear doei0',XX]);
  else
     fid=fopen(['ctd/brutos/doei0',XX,'.cnv'],'r');
     for i = 1:70; fgetl(fid); end
     ctd = fscanf(fid,'%f',[10 inf]);
     ctd =ctd';
  end

  p=ctd(:,1);
min(p)

%    if min(p)>10
%       disp('spike p')
%    end

  T=ctd(:,2);
  S=ctd(:,5);
 
  %%% removendo spikes

  for i=1:length(T)
      if T(i)==-9.990e-29
         T(i)=T(i-1);         
      end
      if S(i)==-9.990e-29
         S(i)=S(i-1);  
      end
  end    

  for i=1:length(T)
      if T(i)==-9.990e-29 
         disp('spike T ')
         disp(num2str(xx(k)))
      end
      if S(i)==-9.990e-29
         disp('spike S ')
         disp(num2str(xx(k)))
      end
  end    

if max(p) > 33
   T = weim(31,'hann',T); T=T';
   S = weim(31,'hann',S); S=S';
else 
   T = weim(7,'hann',T); T=T';
   S = weim(7,'hann',S); S=S';
end
 %    
%    %% separando os dados de descida
%    d = find(p == max(p));  
%    p = p(1:d);
%    T = T(1:d);
%    S = C(1:d);
  
figure(1); 
subplot(121)
plot(T,-p)
subplot(122)
plot(S,-p)
  
disp(num2str(xx(k)))

  
  format long g
  data=[p'; T'; S']; 

  eval(['fid=fopen(''../../../dados/leste1/brutos/lesteI_ctd',XX,'.dat'',''w'');']);

  fprintf(fid,'%12.4f %12.4f %12.4f\n',data);

  clear p T S C CC CCbin CCs binp pbin Tbin Sbin Ss Ts pmax pmin maxp kk
  clear minpl1 lp lx t1 t2 t3 Tgarf Sgarf Siggarf garf Garf garftest
  clear R RR Sbru Tbru pbru Sg TTH Tg XX data minp q sigT sigTH svel
  
end
