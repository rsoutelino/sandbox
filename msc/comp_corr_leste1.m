%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CALCULO DA FUNCAO DE CORRELACAO PARA 
%    OS DADOS HIDROGRAFICOS DA OCEANO LESTE II
%               Rafael Soutelino
%            Mestrado - IOUSP - 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc

% leitura dos dados	

gpan=[];lon=[];lat=[];
for k=1:12
   eval(['load ../mat/lesteI_agp_rad',num2str(k),'.mat']);
   gpan = [gpan agp(1,:)];
   lon = [lon lons];
   lat = [lat lats];
end
 
f0=sw_f(-15);
psig=gpan/f0; psig=-psig';
psig=psig-mean(psig);
clear agp gpan k f0

% montando vetores com as distancias entre todas as estacoes
di=[];dj=[];dd=[];
nest=length(lon);

for i=1:nest
   for j=1:nest
      di = [di i];
      dj = [dj j];
      dd = [dd sw_dist([lat(i) lat(j)],[lon(i) lon(j)],'km')];
   end
end
dd=dd'; di=di';dj=dj';

% calculando os coef de correlacao para cada distancia pre-estabelecida

D = 55;

for k=1:floor(max(dd)/D)
   f = find(dd>=D*(k-1) & dd<=D*k);
   vi = di(f);
   vj = dj(f);
   ccorr(k) = corr2(abs(psig(vi)),abs(psig(vj))); 
   xc(k) = D*k; 
end

xc=xc-D/2;


% ajustando a gaussiana

% rebatendo os valores 
xc2 = [-fliplr(xc) xc];
ccorr2 = [fliplr(ccorr) ccorr];

% interpolando para enriquecer a gaussiana
xc3 = xc2(1):1:xc2(end);
ccorr3 = interp1(xc2,ccorr2,xc3,'linear');


% y=c*exp(-(x-mu).^2/(2*sigma.^2))
% set start parameters: c mu sigma
param0 = [1 0 100];

% call fminsearch
par = fminsearch(@gaussfit,param0,optimset('display','on','tolx',1e-6,'tolfun',1e-6),[xc3;ccorr3],0);

% read best parameters	       
c = par(1,1);           
mu = par(1,2);
sigma = par(1,3);

xgauss = xc3;
ygauss = c*exp(-(xgauss-mu).^2/(2*sigma.^2));

disp(['Lc = ',num2str(sqrt(2)*sigma),' km   ou   Lc = ',num2str(sqrt(2)*sigma/1.852/60),' graus']);
disp(['E^2 = ',num2str(1-c)]);

figure(1)
set(gcf,'color','w')
hold on
plot(xc,ccorr,'k*')
plot(xgauss,ygauss,'r','linewidth',2)
plot([0 1200],[0 0],'k--')
axis([0 500 -0.1 1.1])
xlabel('Lag (km)')
ylabel('Correlacao')
grid
set(gca,'plotboxaspectratio',[1 0.5 1])

print -depsc ../figuras/comp_corr_leste1.eps
!epstopdf ../figuras/comp_corr_leste1.eps
!rm -rf ../figuras/comp_corr_leste1.eps





