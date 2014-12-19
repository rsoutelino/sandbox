%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      SECOES VERTICAIS DE ADCP PARA OS DADOS DA OCEANO LESTE II
%       A PARTIR DOS DADOS PROCESSADOS PELO CODAS - EXPERIMENTO 1
%            Rafael Soutelino - Mestrado IOUSP - ago/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% carregando os dados e tomando variaveis de lon,lat,u,v,t
load contour/oeii_full_uv.mat;
load contour/oeii_full_xy.mat;
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lon = xyt(1,:); lon = lon-360; lat = xyt(2,:); t = xyt(3,:); 
mod = sqrt(u.^2+v.^2);

% removendo spikes remanescentes
f = find(mod >= 1.5);
u(f)=nan;v(f)=nan;

% declarando projecao para trabalhar com m_map
m_proj('mercator','long',[-41 -33.7],'lat',[-20.5 -10.2],'on');
load ../../common/etopo2_leste.mat;

% salvando as posicoes que limitam as radiais
rad1 = [122:157];rad2 = [212:-1:174];rad3 = [216:250];rad4 = [291:-1:256];
rad5 = [292:330];rad6 = [392:-1:346];rad7 = [445:481];rad8 = [524:-1:491];
rad9 = [525:554];rad10 = [588:-1:562];rad11 = [589:620];rad12 = [673:-1:629];

% escolhendo a radial
figure;hold on;set(gcf,'color','w')
m_plot(lon,lat,'.')
m_usercoast('../m_files/costa_leste.mat','patch',[0 0 0])
m_text(-35,-11,'Rad 1','fontsize',10,'fontweight','bold')
m_text(-35.5,-12,'Rad 2','fontsize',10,'fontweight','bold')
m_text(-37,-20,'Rad 12','fontsize',10,'fontweight','bold')
m_grid

rad = menu('ESCOLHA A RADIAL','Radial 1','Radial 2','Radial 3','Radial 4','Radial 5','Radial 6','Radial 7','Radial 8','Radial 9','Radial 10','Radial 11','Radial 12')

close(1)

% criando as variaveis para a radial
eval(['u = u(:,rad',num2str(rad),');'])
eval(['v = v(:,rad',num2str(rad),');'])
eval(['x = lon(rad',num2str(rad),');'])
eval(['y = lat(rad',num2str(rad),');'])

% rotacionando os eixos para obter componente normal
[ans,ang] = sw_dist([y(1) y(end)],[x(1) x(end)],'km');
[mod,dir] = uv2intdir(u*0,v,0,ang);
[usec,vsec] = intdir2uv(mod,dir,0,0);

% INTERPOLANDO
% criando grade nova
xi = x(1):0.01:x(end);
zi = z(1):2:z(end);
[xg,zg] = meshgrid(xi,zi);

% redimensionando grade antiga e removendo NaNs
[x2,z2] = meshgrid(x,z);
[l1,l2] = size(x2);
x2 = reshape(x2,1,l1*l2);
z2 = reshape(z2,1,l1*l2);
vsec2 = reshape(vsec,1,l1*l2);
f = find(isnan(vsec2)==0);
x2 = x2(f); z2 = z2(f); vsec2 = vsec2(f);

% interpolando
vi = griddata(x2,z2,vsec2,xg,zg,'linear'); 
vi = smoo2(vi,-9999,10,1);

% plotando
lv = -0.65:0.02:0.4;
load jet2; colormap(jet2);

figure(1)

%  subplot(211)
%  hold on
%  contourf(x,-z,vsec,lv);shading flat; 
%  caxis([lv(1) lv(end)]); cc = colorbar;
%    colormap(jet2);
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1) pos(2)*3.1 pos(3)/2 pos(4)/2.3])
%  axis([x(1) x(end) -300 0])
%  xlabel('Longitude [graus]')
%  ylabel('Profundidade [m]')
%  tit = ['Radial ',num2str(rad),' (Bruta)'];
%  title(tit,'fontweight','bold')
%  pbaspect([1 0.4 1])
%  
%  subplot(212)
hold on;
contourf(xg,-zg,vi,lv);shading flat; 
caxis([lv(1) lv(end)]); cc = colorbar;
  colormap(jet2);
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1) pos(2)*3.1 pos(3)/2 pos(4)/2.3])
axis([x(1) x(end) -300 0])
xlabel('Longitude [graus]')
ylabel('Profundidade [m]')
tit = ['Radial ',num2str(rad)];
title(tit,'fontweight','bold')
pbaspect([1 0.4 1])
%  
%  eval(['print -depsc figuras/sec_adcp_rad',num2str(rad),'.eps'])
%  eval(['!epstopdf figuras/sec_adcp_rad',num2str(rad),'.eps'])






