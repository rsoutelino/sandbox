%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CALCULO DA FUNCAO DE CORRELACAO PARA 
%    OS DADOS HIDROGRAFICOS DA OCEANO LESTE II
%               Rafael Soutelino
%            Mestrado - IOUSP - 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc

% leitura dos dados	
load ../mat/psi_woa_larga.mat
[l1,l2] = size(psigi);
psi = reshape(psigi,1,l1*l2); 

% tirando NaN
f=find(not(isnan(psi)));
psi = psi(f);
xg = xg(f);
yg = yg(f);

% tem q decimar, se nao o loop nao roda...
lon = xg(1:10:end);
lat = yg(1:10:end);
psi = psi(1:10:end);

% montando vetores com as distancias entre todas os pontos
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

D = 70;

for k=1:floor(max(dd)/D)
   f = find(dd>=D*(k-1) & dd<=D*k);
   vi = di(f);
   vj = dj(f);
   ccorr(k) = corr2(abs(psi(vi)),abs(psi(vj))); 
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
axis([0 1000 -0.1 1.1])
title('Funcao de Correlacao Gaussiana - WOA 2001','fontweight','bold')
ylabel('Correlacao')
xlabel('Lag (km)')
t1 = ['Lc = ',num2str(sqrt(2)*sigma/1.852/60),' ^{\circ}'];
t2 = ['E^2 = ',num2str(1-c)];
text(700,0.7,t1,'fontsize',12,'fontweight','bold')
text(700,0.6,t2,'fontsize',12,'fontweight','bold')
grid

print -depsc ../figuras_larga_escala/corr_woa_larga.eps
!epstopdf ../figuras_larga_escala/corr_woa_larga.eps





