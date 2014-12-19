%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PROGRAMA PARA PLOTAR AS DISTRIBUICOES
%   VERTICAIS E CALCULAR TRANSPORTE A PARTIR
%   DA JUNCAO DOS CAMPOS HORIZONTAIS AO
%    USANDO O REFERENCIAMENTO POR ADCP
%               Oceano Leste 2
%         AGO/2007 - Mestrado - IOUSP
%           Rafael Guarino Soutelino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc;warning off
tic

for sec = 1:50

sec

% definindo os limites da area
lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
m_proj('mercator','long',[min(lonlim) max(lonlim)],'lat',[min(latlim) max(latlim)],'on');
load ../common/etopo2_leste.mat;


% carregando a grade
load ../common/seagrid_leste2.dat;
grid = seagrid_leste2;
xg = grid(:,9); 
yg = grid(:,10);

lj = max(grid(:,1));
li = max(grid(:,2));

xgi = reshape(xg,li,lj);
ygi = reshape(yg,li,lj);

clear seagrid_leste2 grid

load ../common/etopo2_leste.mat
load mat/uv_psi_ref_OEII_150m.mat

%  % preparando para plotar, possibilitando a escolha da radial
%  lpsi = -40000:600:30000;
%  inc = ( max(max(psi)) - min(min(psi)) ) / 40;
%  lpsi2 = -70000:inc:70000;
%  int = num2str(100*round(inc/100));
%  fc = 1;
%  
%  figure(6)
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
%  m_contourf(xgi,ygi,psi,lpsi);shading flat;hold on
%  caxis([-30000 30000])
%  %  m_contour(xgi,ygi,psi,lpsi2,'w')
%  m_quiver(xgi,ygi,U*fc,V*fc,0,'k')
%  [c,h] = m_contourf(xb,yb,zb,[-100 0],'k');
%  set(h,'facecolor',[.8 .8 .8])
%  [c,h] = m_contour(xb,yb,zb,[-100 0],'k');
%  set(h,'color',[.8 .8 .8])
%  m_usercoast('costa_leste.mat','patch',[0 0 0])
%  m_quiver(-40,-11,0.5*fc,0,0,'w')
%  text = ['50 cm s^{-1}'];
%  m_text(-40,-11.3,text,'color','w','fontsize',10,'fontweight','bold')
%  %  text = [num2str(pp),' m'];
%  %  m_text(-34.8,-20.1,text,'color','r','fontweight','bold')
%  title('\Psi Geostrofico Referenciado por ADCP - OEII [m^2 s^{-1}]','fontsize',10,'fontweight','bold');
%  text = ['Contornos: ',int,' m^2s^{-1}'];
%  m_text(-40.5,-11.8,text,'color','w','fontsize',8,'fontweight','bold')
%  m_grid
%  %  cc = colorbar;
%  %  pos = get(cc,'Position');
%  for  l = 1:3:length(xgi) 
%      text = ['Secao: ',num2str(l)];
%      m_text(xgi(l,end),ygi(l,end),text,'color','k','fontweight','bold')
%  end
%  %  set(cc,'Position',[pos(1) pos(2)*1.3 pos(3)/2.5 pos(4)/1.1])


%%% COMECA AQUI A SELECAO DA SECAO PARA CALCULO DE TRANSPORTE
%  disp(' ')
%  sec = input('Selecione a secao para efetuar o calculo do transporte: ');
%  disp(' ')

%%% construindo o eixo de distancia da costa:
xsec = xgi(sec,:); ysec = ygi(sec,:);
dist = [sw_dist(ysec,xsec,'km')];
dist = [0 cumsum(dist)];

% encontrando coordenada do ponto da linha de costa alinhado com a secao,
% para atribuir valor zero a ele:
load ../common/costa_leste.mat
xlc = ncst(:,1);
ylc = ncst(:,2);

%  disp(' ')
%  disp('Achando o ponto mais proximo a linha de costa ...........')
%  disp(' ')
%  
%  dd = [];
%  for t = 1:length(xsec)
%      for u = 1:length(xlc)
%          dd(u,t) = sw_dist([ysec(t) ylc(u)],[xsec(t) xlc(u)],'km'); 
%      end
%  end    
%  
%  [i,j] = find(dd == min(min(dd)));
%  dist = dist - dist(j)*ones(size(dist));


%%% criando topografia para a secao a partir do etopo2

latbt = mean(ysec);
blat = yb(1,:);
b = near(blat,latbt,1);
xbt = xb(:,b);
ybt = yb(:,b);
zbt = zb(:,b);

zbt = interp1(xbt,zbt,[xsec(1):0.01:xsec(end)],'cubic');
ybt = interp1(xbt,ybt,[xsec(1):0.01:xsec(end)],'linear');
xbt = [xsec(1):0.01:xsec(end)];

%%% construindo a secao a partir das distribuicoes horizontais
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Construindo secao ... vai demorar....')
disp(' ')

vsec = [];
cont = 1;

inc = 20; % resolucao vertical
for n = 10:inc:2500
	nn=num2str(n);
        eval(['load mat/uv_psi_ref_OEII_',nn,'m.mat'])
        vsec(cont,:) = V(sec,:); 
        cont = cont + 1;
end

% alisando para retirar descontinuidades
%  vsec = smoo2(vsec,-9999,5,2);

% construindo eixo vertical:
prof = 1:inc:n-1;


figure(2)
set(gcf,'Color',[1 1 1]);
lv = [-0.5:0.005:0.5];
contourf(xsec,-prof,vsec,lv); shading flat; hold on; caxis([lv(1) lv(end)]);
%  if max(max(abs(vsec))) <= 0.2
%     [c,h]=contour(xsec,-prof,vsec,[-0.8:0.1:0.8],'k');
%  else 
%     [c,h]=contour(xsec,-prof,vsec,[-0.8:0.2:0.8],'k');
%  end
axis([min(min(xsec)) max(max(xsec)) -1000 0])
%  clabel(c,h,'labelspacing',500);
set(gca,'plotboxaspectratio',[1 .7 1])
fill([min(xbt) xbt],[min(zbt) zbt],[0 0 0])
latsec = -(round(ysec(1)*10))./10;
tit = ['Vel.Geo.Ref. ADCP - ',num2str(latsec),'^\circS [m s^{-1}]'];
t = title(tit,'fontweight','bold','fontsize',12);
xlabel('[^\circ W]')
ylabel('Profundidade [m]')
cc = colorbar;
pos = get(cc,'Position');
set(cc,'Position',[pos(1) pos(2)*1.9 pos(3)/2.5 pos(4)/1.33])
load /home/rafaelgs/mestrado/proc/adcp/exper1/jet2
colormap(jet2)
drawnow


   eval(['print -depsc figuras/sec_OEII_AO_',num2str(latsec),'S.eps'])
   eval(['print -dpng figuras/sec_OEII_AO_',num2str(latsec),'S.png'])
   eval(['!epstopdf figuras/sec_OEII_AO_',num2str(latsec),'S.eps'])

clear
close all

end
