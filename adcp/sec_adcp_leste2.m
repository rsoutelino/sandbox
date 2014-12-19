%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    SECOES VERTICAIS DE ADCP PARA OS DADOS DO CRUZEIRO PRO-ABROLHOS I
%       A PARTIR DOS DADOS PROCESSADOS PELO CODAS - EXPERIMENTO 1
%            Rafael Soutelino - Mestrado IOUSP - ago/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off

% configuracoes---------------------------------
alis = 's';
cruz = 'oeii'
CRUZ = 'OEII'
lonlim = [-41 -33.7];
latlim = [-20.5 -10.2];
maperr = 's';
cc = 'n';
transp = 'n';
%-----------------------------------------------

% carregando os dados e tomando variaveis de lon,lat,u,v,t
eval(['load ../exper1/contour/',cruz,'_uv.mat']);
eval(['load ../exper1/contour/',cruz,'_xy.mat']);
u = uv(:,1:2:end-1); v = uv(:,2:2:end); 
lon = xyt(1,:); lon = lon-360; lat = xyt(2,:); t = xyt(3,:); 
mod = sqrt(u.^2+v.^2);

% removendo spikes remanescentes
f = find(mod >= 1.5);
u(f)=nan;v(f)=nan;

% declarando projecao para trabalhar com m_map
m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');
load ../../common/etopo2_leste.mat;

% salvando as posicoes que limitam as radiais
rad1 = [213:278]; rad2 = [371:-1:303]; rad3 = [379:442]; rad4 = [526:-1:457];
rad5 = [528:607]; rad6 = [725:-1:628]; rad7 = [827:894];
rad8 = [977:-1:915]; rad9 = [978:1039]; rad10 = [1112:-1:1054]; 
rad11 = [1113:1176]; rad12 = [1295:-1:1196];
% escolhendo a radial
figure;hold on;set(gcf,'color','w')
p = m_plot(lon,lat,'.');
set(p,'color',[.5 .5 .5])
m_usercoast('../../common/costa_leste.mat','patch',[0.542 0.422 0.000])
m_text(-34.9535,-10.9300,'1','fontsize',10,'fontweight','bold')
m_text(-35.8242,-11.7799,'2','fontsize',10,'fontweight','bold')
m_text(-36.4875,-12.5677,'3','fontsize',10,'fontweight','bold')
m_text(-37.2545,-13.3347,'4','fontsize',10,'fontweight','bold')
m_text(-37.1509,-14.3090,'5','fontsize',10,'fontweight','bold')
m_text(-37.0680,-15.3248,'6','fontsize',10,'fontweight','bold')
m_text(-36.8814,-16.0296,'7','fontsize',10,'fontweight','bold')
m_text(-37.0058,-16.7759,'8','fontsize',10,'fontweight','bold')
m_text(-36.6741,-17.4186,'9','fontsize',10,'fontweight','bold')
m_text(-36.4875,-17.9783,'10','fontsize',10,'fontweight','bold')
m_text(-36.6741,-18.6831,'11','fontsize',10,'fontweight','bold')
m_text(-37.3997,-19.5123,'12','fontsize',10,'fontweight','bold')
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
[mod,dir] = uv2intdir(u,v,0,ang);
[usec,vsec] = intdir2uv(mod,dir,0,0);

%  if rad==1 | rad==2 | rad==3 | rad==6;
%     vsec = vsec;
%  elseif rad==4
%     vsec = -usec;
%  else
%     vsec = usec;
%  end

% criando eixo de distancia along track
dist = sw_dist(x,y,'km');
x = [0 cumsum(dist)];

% escolhendo prof max da secao 
imagesc(x,z,vsec);colorbar('horiz');hold on
title('Clique na regiao de alcance maximo do ADCP')
f = ginput(1);
pmax = round(f(2))

%  if isempty(f)==0
%     f = f(:,1);
%     f = round(f);
%     vsec(:,f) = nan;
%  end
close(1)

% criando mascara para os dados
[a,b] = size(vsec);

zb = [];
for k = 1:b;
    perf = vsec(:,k);
    f1 = find(isnan(perf) == 1);
    if length(f1) == a
       zb(k) = pmax;
    else
       f2 = [0 ; diff(f1)];
       f3 = find(f2 == 1);
       zb(k) = z(f1(f3(1)));
    end   
end

% INTERPOLACAO POR ANALISE OBJETIVA

% criando grade nova
xi = min(x):2:max(x);
zi = min(z):2:max(z);
[xg,zg] = meshgrid(xi,zi);

% alisando mascara de topografia
zbi = interp1(x,zb,xi);
zbi = weim(11,'hann',zbi);

% redimensionando grade antiga e removendo NaNs
[x2,z2] = meshgrid(x,z);
[l1,l2] = size(x2);
x2 = reshape(x2,1,l1*l2);
z2 = reshape(z2,1,l1*l2);
vsec2 = reshape(vsec,1,l1*l2);
f = find(isnan(vsec2)==0);
x2 = x2(f); z2 = z2(f); vsec2 = vsec2(f);

% redimensionando grade nova para entrar na objmap.m
[L1,L2] = size(xg);
xg2 = reshape(xg,1,L1*L2);
zg2 = reshape(zg,1,L1*L2);

% escolhendo parametros de interpolacao
lx = 50;
lz = 50; 
E = 0.02;

% implementando condicoes de contorno
% usando topografia oriunda das estacoes CTD

%  eval(['load mat/topo_rad',num2str(rad),'.mat;']);
%  
%  if cc == 's'
%     bx2 = bx(find(bz >= -500));
%     bz2 = bz(find(bz >= -500));
%     x2 = [x2 bx2/1000];
%     z2 = [z2 -bz2];
%     vsec2 = [vsec2 zeros(size(bx2))];
%  end

% interpolando
[vo,er] = objmap(x2',z2',vsec2',xg2',zg2',[lx lz],E); 
vo = reshape(vo,L1,L2);
er = reshape(er,L1,L2);
er = 100*sqrt(er);

% plotando
lv1 = -0.8:0.005:0.5;
lv2 = -1:0.1:1;
load jet2;


figure(1);
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

subplot(211)
contourf(xg,-zg,vo,lv1);shading flat; 
hold on
[c,h] = contour(xg,-zg,vo,lv2,'k');
clabel(c,h);
fill([min(xi) xi max(xi)],-[max(zbi)+100 zbi max(zbi)+100],[1 1 1])
%  fill(bx/1000,bz,[0 0 0])
caxis([lv1(1) lv1(end)]); cc = colorbar;
  colormap(jet2);
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1) pos(2) pos(3)/2 pos(4)])
axis([min(min(xg)) max(max(xg)) -300 0])
xlabel('Distancia ao longo da Radial [km]')
ylabel('Profundidade [m]')
tit = ['Radial ',num2str(rad),' - Velocidades Normais'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])

zb2 = zb;
load ../../common/etopo2_leste.mat;

subplot(212)
[c,h] = m_contour(xb,yb,zb,[-200 -1000],'k'); hold on
clabel(c,h,'labelspacing',500);
p = m_plot(lon,lat,'.k','markersize',6);
set(p,'color',[.5 .5 .5])
eval(['m_plot(lon(rad',num2str(rad),'),lat(rad',num2str(rad),'),''.r'',''markersize'',6)'])
m_usercoast('../../common/costa_leste.mat','patch',[0.542 0.422 0.000])
m_text(-34.9535,-10.9300,'1','fontsize',10,'fontweight','bold')
m_text(-35.8242,-11.7799,'2','fontsize',10,'fontweight','bold')
m_text(-36.4875,-12.5677,'3','fontsize',10,'fontweight','bold')
m_text(-37.2545,-13.3347,'4','fontsize',10,'fontweight','bold')
m_text(-37.1509,-14.3090,'5','fontsize',10,'fontweight','bold')
m_text(-37.0680,-15.3248,'6','fontsize',10,'fontweight','bold')
m_text(-36.8814,-16.0296,'7','fontsize',10,'fontweight','bold')
m_text(-37.0058,-16.7759,'8','fontsize',10,'fontweight','bold')
m_text(-36.6741,-17.4186,'9','fontsize',10,'fontweight','bold')
m_text(-36.4875,-17.9783,'10','fontsize',10,'fontweight','bold')
m_text(-36.6741,-18.6831,'11','fontsize',10,'fontweight','bold')
m_text(-37.3997,-19.5123,'12','fontsize',10,'fontweight','bold')
m_grid

eval(['print -depsc ../figuras/sec_adcp_rad',num2str(rad),'.eps'])
eval(['!epstopdf ../figuras/sec_adcp_rad',num2str(rad),'.eps'])

if maperr == 's'
figure(2);
set(gcf,'Color',[1 1 1])
contourf(xg,-zg,er,[1:0.1:30]);shading flat; 
hold on
%  [c,h] = contour(xg,-zg,er,[0:1:20],'k');
%  clabel(c,h);
fill([min(xi) xi max(xi)],-[max(zbi)+100 zbi max(zbi)+100],[1 1 1])
%  fill(bx/1000,bz,[0 0 0])
caxis([1 30]); cc = colorbar;
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1) pos(2)*2.3 pos(3)/2 pos(4)/1.5])
axis([min(min(xg)) max(max(xg)) -300 0])
xlabel('Distancia ao longo da radial [km]')
ylabel('Depth [m]')
tit = ['Radial ',num2str(rad),' - Erro de Interpolacao'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])

eval(['print -depsc ../figuras/err_sec_adcp_rad',num2str(rad),'.eps'])
eval(['!epstopdf ../figuras/err_sec_adcp_rad',num2str(rad),'.eps'])

end

% plotando linear pra conferir-------------------------------------------
% figure(3);
%  set(gcf,'color','w');
%  contourf(x,-z,vsec,lv1);shading flat; 
%  hold on
%  [c,h] = contour(x,-z,vsec,lv2,'k');
%  caxis([lv1(1) lv1(end)]); cc = colorbar;
%  fill([min(x) x max(x)],-[max(zb2)+100 zb2 max(zb2)+100],[1 1 1])
%    colormap(jet2);
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1) pos(2)*2.4 pos(3)/2 pos(4)/1.6])
%  axis([x(1) x(end) -500 0])
%  xlabel('Distance along ship track [km]')
%  ylabel('Depth [m]')
%  tit = ['Transect ',num2str(rad)];
%  title(tit,'fontweight','bold')
%  pbaspect([1 0.6 1])

if transp == 's'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                VOLUME TRANSPORT CALCULATION                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULA OS VALORES DE VELOCIDADE NO CENTRO DA CELULA DE GRADE %%%
xg2 = xg;

vom=0.5*(vo(:,2:end)+vo(:,1:end-1));
%  vom=0.5*(vo(2:end,:)+vo(1:end-1,:));

%%% MATRIZ DE dx's E dz's %%%

dx=xg2(:,2:end)-xg2(:,1:end-1); 
dz=zg(2:end,:)-zg(1:end-1,:); dz=-dz;

xm=0.5*(xg2(:,2:end)+xg2(:,1:end-1));
xm=0.5*(xm(2:end,:)+xm(1:end-1,:));
zm=0.5*(zg(:,2:end)+zg(:,1:end-1)); zm = -zm;
zm=0.5*(zm(2:end,:)+zm(1:end-1,:)); 

%%% ELIMINANDO 10 ULTIMOS NIVEIS SIGMA %%%
xm = xm(1:end-0,:);
zm = zm(1:end-0,:);
vom = vom(1:end-1,:);

figure(5)
set(gcf,'color','w')
contourf(xg,-zg,vo,lv1);shading flat; 
hold on
[c,h] = contour(xg,-zg,vo,[0 0],'k');
clabel(c,h);
fill([min(xi) xi max(xi)],-[max(zbi)+100 zbi max(zbi)+100],[1 1 1])
%  fill(bx/1000,bz,[0 0 0])
caxis([lv1(1) lv1(end)]); % cc = colorbar;
  colormap(jet2);
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1) pos(2) pos(3)/2 pos(4)])
axis([min(min(xg)) max(max(xg)) -300 0])
xlabel('Along-track distance [km]')
ylabel('Depth [m]')
tit = ['Transect ',num2str(rad),' - Cross-section velocities'];
title(tit,'fontweight','bold')
pbaspect([1 0.6 1])


zm = -zm;
% selecting the currents
disp(' ')
disp('Selecione os limites laterais e verticais da CB para sul que deseja medir: ')
disp(' ')
sul1 = ginput(4);
fsul1 = find(vom<0 & xm>=sul1(1,1) & xm<=sul1(2,1) & -zm<=sul1(3,2) & -zm>=sul1(4,2));
p = plot(xm(fsul1),-zm(fsul1),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais do 2o fluxo para sul que deseja medir: ')
disp(' ')
sul2 = ginput(4);
fsul2 = find(vom<0 & xm>=sul2(1,1) & xm<=sul2(2,1) & -zm<=sul2(3,2) & -zm>=sul2(4,2));
p = plot(xm(fsul2),-zm(fsul2),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais do 3o fluxo para sul que deseja medir: ')
disp(' ')
sul3 = ginput(4);
fsul3 = find(vom<0 & xm>=sul3(1,1) & xm<=sul3(2,1) & -zm<=sul3(3,2) & -zm>=sul3(4,2));
p = plot(xm(fsul3),-zm(fsul3),'k.','markersize',3)
set(p,'color',[.7 .7 .7])

disp(' ')
disp('Selecione os limites laterais e verticais da CCI para norte que deseja medir: ')
disp(' ')
nor1 = ginput(4);
fnor1 = find(vom>0 & xm>=nor1(1,1) & xm<=nor1(2,1) & -zm<=nor1(3,2) & -zm>=nor1(4,2));
plot(xm(fnor1),-zm(fnor1),'w.','markersize',3)

disp(' ')
disp('Selecione os limites laterais e verticais do 2o fluxo para norte que deseja medir: ')
disp(' ')
nor2 = ginput(4);
fnor2 = find(vom>0 & xm>=nor2(1,1) & xm<=nor2(2,1) & -zm<=nor2(3,2) & -zm>=nor2(4,2));
plot(xm(fnor2),-zm(fnor2),'w.','markersize',3)

disp(' ')
disp('Selecione os limites laterais e verticais do 3o fluxo para norte que deseja medir: ')
disp(' ')
nor3 = ginput(4);
fnor3 = find(vom>0 & xm>=nor3(1,1) & xm<=nor3(2,1) & -zm<=nor3(3,2) & -zm>=nor3(4,2));
plot(xm(fnor3),-zm(fnor3),'w.','markersize',3)

%%% TRANSPORTE DE VOLUME %%%

Tsul1 = sum(vom(fsul1).*dx(fsul1).*dz(fsul1))*1e-6*1e3;
Tsul2 = sum(vom(fsul2).*dx(fsul2).*dz(fsul2))*1e-6*1e3;
Tsul3 = sum(vom(fsul3).*dx(fsul3).*dz(fsul3))*1e-6*1e3;

Tnorte1 = sum(vom(fnor1).*dx(fnor1).*dz(fnor1))*1e-6*1e3;
Tnorte2 = sum(vom(fnor2).*dx(fnor2).*dz(fnor2))*1e-6*1e3;
Tnorte3 = sum(vom(fnor3).*dx(fnor3).*dz(fnor3))*1e-6*1e3;

%%% Criando meios de plotar os valores na propria figura
clear text
format bank

xsul1 = mean(mean(xm(fsul1)));
zsul1 = mean(mean(zm(fsul1)));
Tsul1 = (round(Tsul1*10))/10;
tsul1 = [num2str(abs(Tsul1)),' Sv'];
text(xsul1,-zsul1,tsul1,'fontweight','bold')
vmaxsul1 = round((min(min(vom(fsul1))))*100);
tsul1 = [num2str(vmaxsul1),' cm s^{-1}'];
text(xsul1,-(zsul1+100),tsul1,'fontweight','bold')

xsul2 = mean(mean(xm(fsul2)));
zsul2 = mean(mean(zm(fsul2)));
Tsul2 = (round(Tsul2*10))/10;
tsul2 = [num2str(abs(Tsul2)),' Sv'];
text(xsul2,-zsul2,tsul2,'fontweight','bold')
vmaxsul2 = round((min(min(vom(fsul2))))*100);
tsul2 = [num2str(vmaxsul2),' cm s^{-1}'];
text(xsul2,-(zsul2+100),tsul2,'fontweight','bold')

xsul3 = mean(mean(xm(fsul3)));
zsul3 = mean(mean(zm(fsul3)));
Tsul3 = (round(Tsul3*10))/10;
tsul3 = [num2str(abs(Tsul3)),' Sv'];
text(xsul3,-zsul3,tsul3,'fontweight','bold')
vmaxsul3 = round((min(min(vom(fsul3))))*100);
tsul3 = [num2str(vmaxsul3),' cm s^{-1}'];
text(xsul3,-(zsul3+100),tsul3,'fontweight','bold')

xnorte1 = mean(mean(xm(fnor1)));
znorte1 = mean(mean(zm(fnor1)));
Tnorte1 = (round(Tnorte1*10))/10;
tnorte1 = [num2str(abs(Tnorte1)),' Sv'];
text(xnorte1,-znorte1,tnorte1,'fontweight','bold')
vmaxnorte1 = round((max(max(vom(fnor1))))*100);
tnorte1 = [num2str(vmaxnorte1),' cm s^{-1}'];
text(xnorte1,-(znorte1+100),tnorte1,'fontweight','bold')

xnorte2 = mean(mean(xm(fnor2)));
znorte2 = mean(mean(zm(fnor2)));
Tnorte2= (round(Tnorte2*10))/10;
tnorte2 = [num2str(abs(Tnorte2)),' Sv'];
text(xnorte2,-znorte2,tnorte2,'fontweight','bold')
vmaxnorte2 = round((max(max(vom(fnor2))))*100);
tnorte2 = [num2str(vmaxnorte2),' cm s^{-1}'];
text(xnorte2,-(znorte2+100),tnorte2,'fontweight','bold')

xnorte3 = mean(mean(xm(fnor3)));
znorte3 = mean(mean(zm(fnor3)));
Tnorte3 = (round(Tnorte3*10))/10;
tnorte3 = [num2str(abs(Tnorte3)),' Sv'];
text(xnorte3,-znorte3,tnorte3,'fontweight','bold')
vmaxnorte3 = round((max(max(vom(fnor3))))*100);
tnorte3 = [num2str(vmaxnorte3),' cm s^{-1}'];
text(xnorte3,-(znorte3+100),tnorte3,'fontweight','bold')


print(5,'-dpng',['../figuras/transp_adcp_rad',num2str(rad)])
print(5,'-depsc',['../figuras/transp_adcp_rad',num2str(rad)])

end


