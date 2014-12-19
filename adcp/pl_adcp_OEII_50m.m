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
N = 50; %input('Escolha o nivel para plotar os vetores:   ');
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


f2 = [91:252];
u(f2) = nan; v(f2) = nan;

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
  s=1; % escala para vetores
  pl = 100;
  p = N;
  load ../../common/etopo2_leste.mat;
  hold on;
  m_proj('mercator','long',[-41 -33.7],'lat',[-22.5 -12],'on');
  m_quiver(lon,lat,u*s,v*s,0,'k');
  m_plot(lon(f2),lat(f2),'.w');
  [c,h]=m_contourf(xb,yb,zb,[-pl -pl],'k');
  set(h,'facecolor',[.7 .7 .7]);
  [c,h]=m_contour(xb,yb,zb,[-pl -pl],'k');
  set(h,'color',[.7 .7 .7]);
  m_usercoast('costa_leste.mat','patch',[0 0 0])
  m_grid('box','fancy','yaxislocation','left','xaxislocation','bottom');
  text = ['50 cm s^{-1}'];
  m_quiver(-40,-13,0.5*s,0,0,'w')
  m_text(-40,-13.3,text,'color','w','fontsize',10,'fontweight','bold')
  text = [num2str(p),' m'];
  m_text(-34.8,-22,text,'color','k','fontweight','bold')
  title('ADCP Velocities - OEII','fontsize',10,'fontweight','bold');
  hold on
   
  
%  eval(['print -depsc ../figuras/adcp_',num2str(N),'m.eps'])
%  eval(['!epstopdf ../figuras/adcp_',num2str(N),'m.eps'])

  


