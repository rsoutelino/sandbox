%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	PLOTA OS PONTOS DE AMOSTRAGEM 
%     DE CORRENTOGRAFO, MAREGRAFO, CTD E
%         ESTACAO METEOROLOGICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all; clc

format bank

load batimetria/estacoes.dat;
xest=estacoes(:,1);
yest=estacoes(:,2);

lonfund = -(45 + 24.6/60);
latfund = -(23 + 49.6/60);

lonmare = -(45 + 25.2/60);
latmare = -(23 + 49.6/60);

lonmet = - 45.4223333; 
latmet = -23.826333333; 

load batimetria/batimetria_GMS.dat;
ssb = batimetria_GMS;
lon=ssb(:,1);
lat=ssb(:,2);
z=ssb(:,3);

lonlim=[min(lon) max(lon)];
latlim=[min(lat) max(lat)];

m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');

%%% CRIANDO GRADE PARA INTERPOLAR %%%%%%%%

x=min(min(lon)):0.0003:max(max(lon));
y=min(min(lat)):0.0003:max(max(lat));

[xc,yc]=meshgrid(x,y);
[x,y,zc]=griddata(lon,lat,z,xc,yc); % interpolando

[xc,yc]=m_ll2xy(xc,yc);



figure(2)
  set(gcf,...
        'Color',[1 1 1],... 
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

  contourf(xc,yc,zc,min(min(zc)):max(max(zc)));shading flat; 
  
  caxis([-100 200])
  cc = colorbar;
  pos = get(cc,'Position');
  set(cc,'Position',[pos(1)/1.2 pos(2) pos(3)/2.5 pos(4)])
  hold on
%    colorbar('horiz')
  m_usercoast('ctd/costa_ssebastiao.mat','patch',[.6 .6 .6],'LineStyle','-');
%    m_gshhs_f('patch',[.6 .6 .6],'LineStyle','-');
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','right','xaxislocation','top','fontsize',10);
  m_grid('box','fancy','xtick',6,'ytick',6,'yaxislocation','left','xaxislocation','bottom','fontsize',10);
  m_plot(xest,yest,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',8,'MarkerFaceColor',...
       'k','MarkerEdgeColor','w');
%    m_plot(lonfund,latfund,'Color','r','LineStyle','none','LineWidth',...
%         1,'marker','o','markersize',8,'MarkerFaceColor',...
%         'r','MarkerEdgeColor','w');

%  lonmet = - 45.4223333; latmet = -23.826333333;   

%    m_plot(lonmet,latmet,'Color','r','LineStyle','none','LineWidth',...
%         1,'marker','o','markersize',8,'MarkerFaceColor',...
%         'k','MarkerEdgeColor','w');
%    m_plot(lonmare,latmare,'Color','r','LineStyle','none','LineWidth',...
%         1,'marker','o','markersize',8,'MarkerFaceColor',...
%         'y','MarkerEdgeColor','k');
  m_plot(-45.5096970774099,-23.7166244755767,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',8,'MarkerFaceColor',...
       'k','MarkerEdgeColor','w');
  m_text(-45.5026427671442,-23.7166244755767,'Estacao Meteorologica','fontsize',12,'fontweight','bold','color','y');
  m_plot(-45.5096970774099,-23.7281813598226,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',8,'MarkerFaceColor',...
       'r','MarkerEdgeColor','w');
  m_text(-45.5026427671442,-23.7298808152191,'Fundeio - Correntografo','fontsize',12,'fontweight','bold','color','y');
  m_plot(-45.5096970774099,-23.7421162398812,'Color','r','LineStyle','none','LineWidth',...
       1,'marker','o','markersize',8,'MarkerFaceColor',...
       'y','MarkerEdgeColor','k');
  m_text(-45.5026427671442,-23.7438155135424,'Maregrafo','fontsize',12,'fontweight','bold','color','y');

[x,y]=ginput(5);
[x,y]=m_xy2ll(x,y);

m_text(x(1),y(1),'1','fontsize',12,'fontweight','bold','color','y');
m_text(x(2),y(2),'2','fontsize',12,'fontweight','bold','color','y');
m_text(x(3),y(3),'3','fontsize',12,'fontweight','bold','color','y');
m_text(x(4),y(4),'4','fontsize',12,'fontweight','bold','color','y');
m_text(x(5),y(5),'5','fontsize',12,'fontweight','bold','color','y');


%   print -depsc grade.eps
%  !epstopdf grade.eps
stop

format long g

[x0,y0] = ginput(1);
plot(x0,y0,'k.','markersize',10)
 
inc=-0.0005;

xeixo1=[x0 x0]; yeixo1=[y0+inc y0-inc];
xeixo2=[x0+inc x0-inc]; yeixo2=[y0 y0];
plot(xeixo1,yeixo1,'--k')
plot(xeixo2,yeixo2,'--k')

% criando o eixo rotacionado
inc2=-0.0004

xrot1=[x0-inc2 x0+inc2]; yrot1=[y0+inc2 y0-inc2];
plot(xrot1,yrot1,'y','linewidth',2)
xrot2=[x0-inc2 x0+inc2]; yrot1=[y0 y0];
xrot2=[x0-inc2 x0+inc2]; yrot2=[y0-inc2 y0+inc2];
plot(xrot2,yrot2,'y','linewidth',2)

g=ginput(1);
text(g(1),g(2),'E')
g=ginput(1);
text(g(1),g(2),'N')
g=ginput(1);
text(g(1),g(2),'cross-channel','color','y','fontweight','bold')
g=ginput(1);
text(g(1),g(2),'along-channel','color','y','fontweight','bold')
g=ginput(1)
text(g(1),g(2),'45^\circ','color','y','fontsize',8)


