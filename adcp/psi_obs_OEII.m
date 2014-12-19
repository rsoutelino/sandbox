%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FUNCAO DE CORRENTE OBSERVADA PARA OS DADOS DA OCEANO LESTE II
%       A PARTIR DOS DADOS PROCESSADOS PELO CODAS - EXPERIMENTO 1
%            Rafael Soutelino - Mestrado IOUSP - ago/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% carregando os dados e tomando variaveis de lon,lat,u,v,t
load contour/oeii_full_uv.mat;
load contour/oeii_full_xy.mat;
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lon = xyt(1,:); lon = lon-360; lat = xyt(2,:); t = xyt(3,:);  

% declarando projecao para trabalhar com m_map
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../../common/etopo2_leste.mat;

% escolhendo profundidade a ser plotada
p = input('Escolha a prof de interesse (20 a 460 de 10 em 10): ')
n = find(zc==p);
u = u(n,:); v = v(n,:);
mod = sqrt(u.^2+v.^2);

% removendo NaNs e spikes remanescentes
f = find(mod >= -1.5 & mod <= 1.5 & isnan(mod)==0);
u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);

% preparando para a analise objetiva vetorial
load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

l1 = max(grid(:,1));
l2 = max(grid(:,2));

xgi = reshape(xg,l2,l1);
ygi = reshape(yg,l2,l1);

lc = 1; % em graus, entao, precisa passar as velocidades para grau/s
E = 0.02;

u = u*8.9993e-06; % fator de conversao para grau/s (9e-6)
v = v*8.9993e-06;

% analise objetiva vetorial para obter psi observado
[psiob] = vectoa(xg,yg,lon,lat,u,v,lc,E,0);

% aplicando condicao de contorno de dirichilet
figure(1)
if p <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-p -p]);
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
%  xcont = weim(31,'hann',xcont);
%  ycont = weim(31,'hann',ycont);
xcont=xcont';
ycont=ycont';

close(1); 

xf=[min(xcont) xcont];
yf=[max(ycont) ycont];
[xf,yf] = m_ll2xy(xf,yf,'patch');

psiob = psiob*(111120^2);
psiob = psiob-mean(psiob);

%% aplicando condicao de contorno em psiob
xg2 = [xg' xcont]; yg2 = [yg' ycont];
psiob = [psiob' zeros(size(xcont))];
%  xg2 = xg; yg2 = yg;

[psiob] = scaloa(xg',yg',xg2,yg2,psiob,lc,E);
psiob = reshape(psiob,l2,l1);
[uo,vo] = psi2uv(xgi,ygi,psiob);

subplot(121)
m_contourf(xgi,ygi,psiob,20);shading flat;hold on
m_contour(xgi,ygi,psiob,20,'w'); 
%  m_quiver(xgi,ygi,uo,vo,3,'k');
m_quiver(lon,lat,u,v,2,'k'); 
fill(xf,yf,[.8 .8 .8])
m_usercoast('../m_files/costa_leste.mat','patch',[0 0 0])
m_grid

subplot(122)
m_contourf(xgi,ygi,psiob,20);shading flat;hold on
m_contour(xgi,ygi,psiob,20,'w'); 
m_quiver(xgi,ygi,uo,vo,3,'k');
%  m_quiver(lon,lat,u,v,2,'k'); 
fill(xf,yf,[.8 .8 .8])
m_usercoast('../m_files/costa_leste.mat','patch',[0 0 0])
m_grid

%  print -depsc teste.eps

Uadcp = uo;
Vadcp = vo;

save ../../adcp/mat/adcp_codas.mat  Uadcp Vadcp













