%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotando os dados do ADCP processados pelo CODAS %%
%%         Rafael Soutelino - set/2007              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;warning off

% configuracoes---------------------------------
alis = 's';
jan = 11;
ans = 1;
cruz = 'oeii'
sec = 'n'
%-----------------------------------------------

% carregando os dados
eval(['load ../exper1/contour/',cruz,'_uv.mat;']);
eval(['load ../exper1/contour/',cruz,'_xy.mat;']);

u = uv(:,1:2:end-1);
v = uv(:,2:2:end);
lon = xyt(1,:)-360; lat = xyt(2,:);
t = xyt(3,:);

%% eliminando valores de magnitude maior que 1m/s
mod = sqrt(u.^2+v.^2);
f = find(mod > 1);
u(f)=nan; v(f)=nan; mod(f)=nan;

% plotando vetores em uma camada escolhida
disp('Dados de 20m a 372m, de 8/8m)');
N = input('Escolha o nivel para plotar os vetores:   ');
n = near(zc,N,1);

if alis == 's'
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

  set(gcf,'Color','w');axes('color','none'); hold on;
  lonlim = [-41 -33.7]; latlim = [-20.5 -10.2];
  m_proj('mercator','long',[lonlim(1) lonlim(2)],'lat',[latlim(1) latlim(2)],'on');
  [data,lon_bat,lat_bat]=m_tbase([lonlim(1) lonlim(2) latlim(1) latlim(2)]);
  data2=smoo2(data,-9999,9,3.5);
  bat=[-100 -300 -1000 -4000];
  load /usr/local/matlab/toolbox/lado/m_map/ocean_colormap
  m_contourf(lon_bat,lat_bat,data2,50);shading flat;hold on;colormap(ocean_colormap)
  [c,h]=m_contour(lon_bat,lat_bat,data2,bat,'k'); 
  clabel(c,h,'VerticalAlignment','middle','fontsize',10,'labelspacing',400);
  m_usercoast('costa_leste.mat','patch',[0 0.2 0],'LineStyle','-');
  m_grid('box','fancy')
  s=1; % escala para vetores
  m_text(-40,-11.5,'50 cm/s','Fontsize',10,'FontWeight','Bold','Color','y');
  m_quiver(-40,-11.2,0.5*s,0*s,0,'y')
  tex = [num2str(N),' m'];
  m_text(-34.8,-20,tex,'color','r','fontweight','bold')
  hold on
  m_quiver(lon,lat,u*s,v*s,0,'y'); 
  stop
eval(['print -depsc ../figuras/adcp_',num2str(N),'m.eps'])
eval(['!epstopdf ../figuras/adcp_',num2str(N),'m.eps'])

  


