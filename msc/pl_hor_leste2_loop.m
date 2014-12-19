%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PROGRAMA PARA PLOTAR AS DISTRIBUICOES
%        HORIZONTAIS - Oceano Leste 2
%         out/2006 - Mestrado - IOUSP
%           Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Plota distribuições horizontais de T e S, 
%    através de uma grade curvilinear, usando 
%    interpolacao linear e interpolacao 
%    objetiva
% - mostra a grade e sua ortogonalidade
% calcula funcao de corrente geostrofica

    
clear all;close all;clc;warning off

cc = 'n';

for n = 10:10:2500

% definindo os limites da area
lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
%  m_gshhs_f('save','costa_leste.mat')

% lendo as posicoes das estacoes

load posicoes_leste2.dat; pb=posicoes_leste2;
latg=pb(:,2); latm=pb(:,3); latctd=latg+latm/60; latctd=-latctd; latctd=latctd'; 
long=pb(:,4); lonm=pb(:,5); lonctd=long+lonm/60; lonctd=-lonctd; lonctd=lonctd'; 

clear latg latm long lonm posicoes_leste2

%%% carregando os dados hidrograficos %%%%%

%  n=input('Escolha a profundidade: '); 


nn=num2str(n)
% carregando matrizes com os dados hidrograficos
gpan=[];T=[];S=[];lon=[];lat=[];
for k=1:12
   eval(['load ../mat/lesteII_agp_rad',num2str(k),'.mat']);
   eval(['load ../mat/lesteII_T_rad',num2str(k),'.mat']);
   eval(['load ../mat/lesteII_S_rad',num2str(k),'.mat']);

   gpan = [gpan agp(n,:)];
   lon = [lon lons];
   lat = [lat lats];  
   T = [T Ti(n,:)];
   S = [S Si(n,:)];
end

clear agp Ti Si lons lats

%%% calculando PSI
f0=sw_f(-15);
psig=gpan/f0; psig=-psig;
psig=psig-mean(psig);
psigctd = psig;

%%% carregando a grade construida no seagrid

load ../../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

%calculando a ortonogonalidade da grade
[ort,xo,yo]=orthog(xgi,ygi);

%%% interpolacao linear com grade curvilinear %%%%%%%%%%%%%

disp('Interpolacao Linear.........')
Ti = griddata(lon,lat,T,xgi,ygi);
Si = griddata(lon,lat,S,xgi,ygi);
psigi = griddata(lon,lat,psig,xgi,ygi);
[ui,vi]=psi2uv(xgi,ygi,psigi);

%%% ANALISE OBJETIVA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('Interpolacao Objetiva das variaveis termohalinas')
disp(' ')
xg=xg';yg=yg';

corrlen = 1; 
err = 0.02; 
[To,er] = scaloa(xg,yg,lon,lat,T,corrlen,err); 
[So,er] = scaloa(xg,yg,lon,lat,S,corrlen,err);
Do=sw_dens0(So,To); Do=Do-1000*ones(size(Do));

%%% APLICANDO CONDICOES DE CONTORNO PARA PSI ********************************
%%% DIRICHILET

%  load ../../levitus/etopo2_E_751x901.xyz;
%  x=etopo2_E_751x901(:,1);
%  y=etopo2_E_751x901(:,2);
%  z=etopo2_E_751x901(:,3);
%  clear etopo2_E_751x901
%  x = reshape(x,751,901);
%  y = reshape(y,751,901);
%  z = reshape(z,751,901);
%  fx = find(x(:,1) >= -43  & x(:,1) <= -30);
%  fy = find(y(1,:) >= -23 & y(1,:) <= -10);
%  xb = x(fx,fy);
%  yb = y(fx,fy);
%  zb = z(fx,fy);
%  save etopo2_leste2.mat xb yb zb

load ../../common/etopo2_leste.mat;

% buscando valor de lon e lat das isobatas de ?? m
figure(1)
if n <= 100
   [c,h] = contour(xb,yb,zb,[-100 -100]);
else 
   [c,h] = contour(xb,yb,zb,[-n -n]);
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

% enxertando valores da interpolacao linear
xadd=[]; yadd=[]; psiadd=[];
for k = 40:45
   xadd = [xadd xgi(k,:)];
   yadd = [yadd ygi(k,:)];
   psiadd = [psiadd psigi(k,:)];
end

% enxertando valores da isobata nos vetores
if cc == 's'
   lon = [lon xcont xadd]; lat = [lat ycont yadd];
   psig = [psig zeros(size(xcont)) psiadd];
end
% retirando eventuais NaNs

psig = psig(find(not(isnan(psig))));
lon = lon(find(not(isnan(psig))));
lat = lat(find(not(isnan(psig))));

disp('Interpolacao Objetiva de Psi')
[psigo,er] = scaloa(xg,yg,lon,lat,psig,corrlen,err);
er=100*sqrt(er); 

To=reshape(To,li,lj);
So=reshape(So,li,lj);
Do=reshape(Do,li,lj);
psigo=reshape(psigo,li,lj);
er=reshape(er,li,lj);

% calculando componentes u e v
[uo,vo]=psi2uv(xgi,ygi,psigo);


%%% PLOTANDO OS RESULTADOS **********************************************************
isob=[-100 -100];
% Distribuicoes termohalinas horizontais:------------------------------------------------------

P = ['ToSoDo'];

lT=26.8:0.05:28.4; %% para superficie %%
lS=36.6:0.03:37.4;  %% para superficie %%
lD=23.3:0.03:24.3; %% para superficie %%
l = ['lTlSlD'];

t1 = ['Temperatura ( \circ C) - OEII - ',nn,' m'];
t2 = ['Salinidade - OEII - ',nn,' m'];
t3 = ['Densidade Potencial (kg m^{-3}) - OEII - ',nn,' m'];
tit = ['t1t2t3'];
%  
%  c=0;
%  for k = 1:3
%    figure(k)
%    set(k,'Color','w');
%    hold on;
%    eval(['[c1,h1] = m_contourf(xgi,ygi,',P(c+k:c+k+1),',',l(c+k:c+k+1),');']); 
%    shading flat; hold on
%    %  clabel(c1,h1); % para aparecer os numeros no mapa %
%  %    [c2,h2] = m_contour(xb,yb,zb,isob,'k');
%  %    clabel(c2,h2,'labelspacing',500)
%    fill(xf,yf,[0.8 0.8 0.8]);
%    m_usercoast('costa_leste.mat','patch',[0 0 0],'LineStyle','-');
%    %  m_grid('box','fancy','yaxislocation','right','xaxislocation','top','fontsize',10);
%    m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%    cc = colorbar;
%    pos = get(cc,'Position');
%    set(cc,'Position',[pos(1)/1.2 pos(2) pos(3)/2.5 pos(4)])
%    hold off;
%    eval(['title(',tit(c+k:c+k+1),',''fontsize'',12,''fontweight'',''bold'')']);
%  drawnow
%       print(k,'-depsc',['../figuras/',P(c+k),'_OEII_',nn,'m']);
%       eval(['!epstopdf ../figuras/',P(c+k),'_OEII_',nn,'m.eps'])
%  c=c+1; 
%  end
%  clear k c

% ***********************************************************************************************
% fator de escala dos vetores
fc = 1; 

% construindo escalas de contorno para psi

lpsi = -40000:800:30000;
inc = ( max(max(psigo)) - min(min(psigo)) ) / 40;
lpsi2 = -30000:inc:30000;
int = num2str(100*round(inc/100));

figure(4)
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

m_contourf(xgi,ygi,psigo,lpsi);hold on;shading flat;
caxis([-30000 30000])
%  m_contour(xgi,ygi,psigo,lpsi2,'w')
q = m_quiver(xgi,ygi,uo*fc,vo*fc,0,'k');
%  set(q,'color',[.3 .3 .3])
fill(xf,yf,[0.8 0.8 0.8]);
%  plot(xf,yf,'color',[.8 .8 .8])
%  m_plot(lonctd,latctd,'w.','markersize',6)
m_usercoast('costa_leste.mat','patch',[0 0 0])
vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
m_quiver(-40,-11,0.2*fc,0,0,'w')
text = ['20 cm s^{-1}'];
m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
text = [nn,' m'];
m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
text = ['Contornos: ',int,' m^2s^{-1}'];
m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',10);
title('\Psi Geostrofico OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
   print(4,'-depsc',['../figuras/psi_OEII_',nn,'m']);
   eval(['!epstopdf ../figuras/psi_OEII_',nn,'m.eps'])
drawnow
%     print(4,'-depsc',['../figuras/psi_OEII_',nn,'m_3gr']);
%     eval(['!epstopdf ../figuras/psi_OEII_',nn,'m_3gr.eps'])

ugo = uo;
vgo = vo;
eval(['save ../mat/uv_psig_AO_OEII_',nn,'m.mat ugo vgo'])

% interpolacao linear	

%  fc = 1; % fator de escala dos vetores
%  %  
%  figure(5)
%  set(gcf,...
%          'Color',[1 1 1],...
%          'InvertHardcopy','on',...
%          'PaperUnits','inches',...
%          'Units','inches',...
%          'PaperOrientation','portrait',...
%          'PaperPosition',[0 0 8.5 11],...
%          'PaperPositionMode','manual',...
%          'PaperType','usletter',...
%          'Position',[.2 .2 8 9],...
%          'ShareColors','off',...
%          'Clipping','on');
%  
%  m_contourf(xgi,ygi,psigi,lpsi);hold on;shading flat;
%  caxis([-30000 30000])
%  m_contour(xgi,ygi,psigi,lpsi2,'w')
%  m_quiver(xgi,ygi,ui*fc,vi*fc,0,'k')
%  fill(xf,yf,[0.8 0.8 0.8]);
%  %  plot(xf,yf,'color',[.8 .8 .8])
%  m_plot(lonctd,latctd,'w.','markersize',6)
%  m_usercoast('costa_leste.mat','patch',[0 0 0])
%  vmax = round(max(max(sqrt((uo*100).^2+(vo*100).^2))));
%  m_quiver(-40,-11,0.2*fc,0,0,'w')
%  text = ['20 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  text = [nn,' m'];
%  m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
%  text = ['Contornos: ',int,' m^2s^{-1}'];
%  m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
%  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom','fontsize',8);
%  title('\Psi Geostrofico OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
%  cc = colorbar;
%  pos = get(cc,'Position');
%  set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])
%  %     print(4,'-depsc',['../figuras/psi_OEII_',nn,'m']);
%  %     eval(['!epstopdf ../figuras/psi_OEII_',nn,'m.eps'])
%  drawnow
%  
%  print -depsc ../../../relatorio1/linear.eps
%  !epstopdf ../../../relatorio1/linear.eps
%  !rm -rf ../../../relatorio1/linear.eps

clear
close all

end

