%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotando os dados do ADCP processados pelo CODAS %%
%%         Rafael Soutelino - jan/2008              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;warning off;

% configuracoes---------------------------------
pathname = ['/home/rafaelgs/doutorado/data/']
alis = 'y';
jan = 11;
ans = 1;
cruz = 'proab1'
CRUZ = 'Proabrolhos 1'
PRINT = 'y';
s = 1; % escala para vetores
%-----------------------------------------------
lonlim = [-45 -30]; latlim = [-24 -9];
m_proj('mercator','long',[lonlim],'lat',[latlim],'on');
m_gshhs_l('save','costa.mat')
[zb,xb,yb] = m_tbase([lonlim latlim]);

eval(['load ',pathname,cruz,'/adcp/',cruz,'_uv.mat;']);
eval(['load ',pathname,cruz,'/adcp/',cruz,'_xy.mat;']);

u = uv(:,1:2:end-1);
v = uv(:,2:2:end);
lon = xyt(1,:)-360; lat = xyt(2,:);
t = xyt(3,:);

%% eliminando valores de magnitude maior que 1m/s
mod = sqrt(u.^2+v.^2);
f = find(mod > 1.5);
u(f)=nan; v(f)=nan; mod(f)=nan;

%% calculando vetores de tempo:
%  year = 2008;
%  [m,d,h,min] = decjul2greg(year,t);
%  
%  if h(end) > 9
%     H = num2str(h(end));
%  else
%     H = ['0,',num2str(h(end))];
%  end
%  
%  if min(end) > 9
%     MIN = num2str(min(end));
%  else
%     MIN = ['0,',num2str(min(end))];
%  end
 

% plotando vetores em uma camada escolhida
disp('Dados de 32 m a 800 m, de 16/16 m)');
N = input('Escolha o nivel para plotar os vetores:   ');
n = near(zc,N,1);

if alis == 'y'
   u = u(n,:); v = v(n,:);
   % removendo NaNs e spikes remanescentes
   f = find(isnan(u)==0);
   u=u(f);v=v(f);lon=lon(f);lat=lat(f);t=t(f);
   u=weim(jan,'hann',u);
   v=weim(jan,'hann',v);
else
   u = u(n,:); v = v(n,:);
end
  

  figure(2)
%  if k == 1
  set(gcf,'Color','w');axes('color','none'); hold on;
  bat1 = [0:100:5000];
  bat2 = [100 200 1000 2500 3000];
  load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap
  cm1 = colormap(ocean_colormap);cm1 = flipud(cm1);
  m_contourf(xb,yb,-zb,bat1); shading flat; colormap(cm1);%colorbar
  [c,h] = m_contour(xb,yb,-zb,bat2,'k');hold on
  c = clabel(c,h,'labelspacing',1000);
  set(h,'color',[.4 .4 .4])
  set(c,'color',[.4 .4 .4],'fontsize',8)
  m_usercoast('costa.mat','patch',[0 0.2 0],'LineStyle','-');
  m_grid('box','fancy')
  m_text(lonlim(1)+1,latlim(1)+1.5,'50 cm/s','Fontsize',10,'FontWeight','Bold','Color','y');
  m_quiver(lonlim(1)+1,latlim(1)+2,0.5*s,0*s,0,'y')
  tex = [num2str(N),' m'];
  m_text(lonlim(1)+1,latlim(2)-1,tex,'color','w','fontweight','bold')
  tit = ['Vetores de Velocidade - ADCP - ',CRUZ];
  title(tit,'fontsize',12,'fontweight','bold')
  hold on
  m_quiver(lon,lat,u*s,v*s,0,'y'); 
if PRINT == 'y'
   print(2,'-depsc',[pathname,'figures/adcpvet_',cruz,'_',num2str(N),'m']);
   eval(['!epstopdf ',pathname,'figures/adcpvet_',cruz,'_',num2str(N),'m.eps'])
end



