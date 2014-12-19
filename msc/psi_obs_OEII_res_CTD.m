%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FUNCAO DE CORRENTE OBSERVADA PARA OS DADOS DO PRO-ABROLHOS I
%       A PARTIR DOS DADOS PROCESSADOS PELO CODAS - EXPERIMENTO 1
%            Rafael Soutelino - Mestrado IOUSP - ago/2007
%     este programa usa os vetores de ADCP na mesma localidade das
%                         estacoes de CTD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% configuracoes ----------------------------------------------------
alis = 'n';
jan = 11;
cruz = 'oeii';
tr = 'oeii';
cc = 's';
maperr = 'n';
%-------------------------------------------------------------------

% carregando os dados e tomando variaveis de lon,lat,u,v,t
eval(['load ../exper1/contour/',cruz,'_full_uv.mat;'])
eval(['load ../exper1/contour/',cruz,'_full_xy.mat;'])
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lon = xyt(1,:); lon = lon-360; lat = xyt(2,:); t = xyt(3,:);  

% declarando projecao para trabalhar com m_map
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../../common/etopo2_leste.mat;

% escolhendo profundidade a ser plotada
p = input('Escolha a prof de interesse (20 a 460 de 10 em 10): ')
n = near(zc,p,1);
u = u(n,:); v = v(n,:);
mod = sqrt(u.^2+v.^2);

% removendo NaNs e spikes remanescentes
f = find(u >= -1.0 & u <= 1.0 & isnan(u)==0);
u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);
f = find(v >= -1.0 & v <= 1.0 & isnan(v)==0);
u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);

if alis == 's'
   u=weim(jan,'hann',u);
   v=weim(jan,'hann',v);
end

% decimando os vetores de ADCP nos pontos mais proximos das
% estacoes CTD:
load ../../hidrografia/posicoes_leste2.dat; pb = posicoes_leste2;
latctd = - (pb(:,2) + pb(:,3)./60);
lonctd = - (pb(:,4) + pb(:,5)./60);

inc=0.05;  
F=[];

for i=1:length(latctd)

    f=find(lat<(latctd(i)+inc) & lat>(latctd(i)-inc) & lon<(lonctd(i)+inc) & lon>(lonctd(i)-inc)); 
    F=[F f];    
end       

u = u(F);
v = v(F);
lat = lat(F);
lon = lon(F);

% acertando o contorno junto aos bancos
if cc == 's'
   u = [u 0 0]; v = [v -0.5 -0.5]; lon = [lon -37.75 -37.25]; lat = [lat -16 -18];
end

% preparando para a analise objetiva vetorial
load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

l1 = max(grid(:,1));
l2 = max(grid(:,2));

xgi = reshape(xg,l2,l1);
ygi = reshape(yg,l2,l1);

ui = griddata(lon,lat,u,xgi,ygi,'linear');
vi = griddata(lon,lat,v,xgi,ygi,'linear');

lc = 1; % em graus, entao, precisa passar as velocidades para grau/s
E = 0.02;

u = u*8.9993e-06; % fator de conversao para grau/s (9e-6)
v = v*8.9993e-06;

% analise objetiva vetorial para obter psi observado
[psiob] = vectoa(xg,yg,lon,lat,u,v,lc,E,0);

if maperr == 's'
  [ans,er] = scaloa(xg,yg,lon,lat,u,lc,E);
  er = 100*sqrt(er);
end

% aplicando condicao de contorno de dirichilet
figure(1)
if p <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
   pl = 100;
else 
   [c,h] = contour(xb,yb,zb,[-p -p]);
   pl = p;
end
xcont = get(h(1),'xdata'); xcont = xcont(find(not(isnan(xcont))));
ycont = get(h(1),'ydata'); ycont = ycont(find(not(isnan(ycont))));
%  xcont = weim(31,'hann',xcont);
%  ycont = weim(31,'hann',ycont);
xcont=xcont';
ycont=ycont';

close(1); 

psiob = psiob*(111120^2);
psiob = psiob-mean(psiob);

%% aplicando condicao de contorno em psiob
if cc == 's'
   xg2 = [xg' xcont]; yg2 = [yg' ycont];
   psiob = [psiob' zeros(size(xcont))];
   [psiob] = scaloa(xg',yg',xg2,yg2,psiob,lc,E);
else
   xg2 = xg'; yg2 = yg';
   psiob = psiob';
end

psiob = reshape(psiob,l2,l1);
[uo,vo] = psi2uv(xgi,ygi,psiob);

if maperr == 's'
   er = reshape(er,l2,l1);
end

s=1;
u = u/8.9993e-06; % reconvertendo para plotar vetores brutos
v = v/8.9993e-06;

%  subplot(121)
%  m_contourf(xgi,ygi,psiob,50);shading flat;hold on
%  m_contour(xgi,ygi,psiob,20,'w'); 
%  %  m_quiver(xgi,ygi,uo,vo,3,'k');
%  m_quiver(lon,lat,u*s,v*s,0,'k'); 
%  [c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
%  set(h,'facecolor',[.9 .9 .9]);
%  [c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
%  set(h,'color',[.9 .9 .9]);
%  m_usercoast('/home/proabrolhos/proabrolhos1/grade/costa_proabro.mat','patch',[0 0 0])
%  m_grid

% construindo escalas de contorno para psi
lpsi = -40000:600:30000;
inc = ( max(max(psiob)) - min(min(psiob)) ) / 40;
lpsi2 = -30000:inc:30000;
int = num2str(100*round(inc/100));

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');

m_contourf(xgi,ygi,psiob,lpsi);hold on; shading flat
caxis([-30000 30000])
%  m_contour(xgi,ygi,psiob,lpsi2,'k'); 
m_quiver(xgi(1:1:end),ygi(1:1:end),uo(1:1:end)*s,vo(1:1:end)*s,0,'k');
[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
m_quiver(-40,-11,0.5*s,0,0,'w')
text = ['50 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [num2str(p),' m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
title('\Psi Observado [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*2.2 pos(3)/2.5 pos(4)/1.5])
   print(1,'-depsc',['../figuras/psiob_OEII_',num2str(p),'m']);
   eval(['!epstopdf ../figuras/psiob_OEII_',num2str(p),'m.eps'])
drawnow

Uadcp = uo; Vadcp = vo;
eval(['save ../mat/psi_obs_',num2str(p),'m.mat Uadcp Vadcp'])

% figura com o mapa de erro da AO
if maperr == 's'
figure(2)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8 9],...
        'ShareColors','off',...
        'Clipping','on');

m_contourf(xgi,ygi,er,50);hold on; shading flat
caxis([5 100])
[c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
set(h,'facecolor',[.7 .7 .7]);
[c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
set(h,'color',[.7 .7 .7]);
m_usercoast('costa_leste.mat','patch',[0 0 0])
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
title('Erro de Interpolacao - ADCP [%]','fontsize',10,'fontweight','bold');
text = [num2str(p),' m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*2.2 pos(3)/2.5 pos(4)/1.5])
   print(2,'-depsc',['../figuras/err_OEII_',num2str(p),'m']);
   eval(['!epstopdf ../figuras/err_OEII_',num2str(p),'m.eps'])
drawnow
end



